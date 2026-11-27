classdef (Sealed) Helmholtz < BaseObject & EventSender
    %HELMHOLTZ Three coils on three supplies, and the field they make at the NVs.
    %
    % ===================================================================
    % THE MODEL -- one line, and everything else follows from it
    % ===================================================================
    %
    %       B = M * I + b0          (B in Gauss, I = [Ix Iy Iz] in A, M 3x3 in G/A)
    %
    % The CURRENTS are the state. B, Brho, Btheta and Bphi are computed from them on every
    % read, so there is nothing to keep in sync and nothing that can disagree.
    %
    % WHICH FRAME B IS IN depends on where M came from:
    %
    %   CALIBRATED (calibCoilMatrix, level 3): M and b0 were fitted to NV spectra, and B is
    %     in the NV CRYSTAL FRAME -- x, y, z along the diamond's cubic axes as the NVs define
    %     them. A field along an NV axis is then exactly theta = 54.7356, phi = 45 (+k*90),
    %     and alignToNV / BtoNV put it there with no scan.
    %   NOMINAL (no calibration loaded): M = diag(magnetic_field_ratio), b0 =
    %     magnetic_field_offset from setupInfo.json, and B is in the coil frame as those
    %     numbers describe it -- with no cross-talk and no measured scale. Only roughly the
    %     crystal frame. isCalibrated says which; the GUI shows it.
    %
    % ===================================================================
    % ANGLES -- the physics (ISO 80000-2) convention
    % ===================================================================
    %
    %   Btheta   POLAR angle, measured from +z, 0..180 deg
    %   Bphi     AZIMUTH, from +x toward +y, -180..180 deg
    %
    %   B = Brho * [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    %
    % NB the class used to name these the other way round (theta azimuthal, phi from z).
    % Anything written against the old naming that reads or sets an angle by name needs
    % the two swapped. Saved states held only the cartesian B, so they are unaffected.
    %
    % ===================================================================
    % LIMITS
    % ===================================================================
    %
    % Each supply runs from minCurrent to maxCurrent (0 to 5.2 A on setup 7 -- unipolar,
    % so no coil can reverse its field). The fields the coils can make are therefore the
    % parallelepiped {M*I + b0 : I in that box}, NOT a box in B, and certainly not the
    % same box for every direction. limits(what) answers "how far can this one coordinate
    % go with the other two held where they are", exactly:
    %
    %   Bx, By, Bz, Brho   the currents are LINEAR in the coordinate, so each supply limit
    %                      is one bound on it and the answer is their intersection
    %   Btheta, Bphi       each current is a*cos(t) + c*sin(t) + e, so each supply limit is
    %                      crossed at closed-form angles; the answer is the stretch between
    %                      crossings that contains the present angle
    %
    % ===================================================================
    % STATE ON DISK -- expProperties.json, key "Helmholtz"
    % ===================================================================
    %
    % h.save(jsonFile) writes the currents, the field they make, the step sizes and the
    % coil calibration into one top-level "Helmholtz" block of expProperties.json, through
    % expPropsWrite (backup first, verified after). Experiments ignore it: loadJsonParams
    % only copies keys that match one of the experiment's own properties, and none is called
    % Helmholtz. h.load(jsonFile) reads it back. An old helmholtzState.mat still loads: its
    % B was in the nominal coil frame, so it is converted to the currents that made it.

    properties (Constant)
        AXIS = ['X', 'Y', 'Z'];
        CONTROL_STATE = {'ON', 'AUTO', 'OFF'};
        NAME = 'Helmholtz';
        JSON_KEY = 'Helmholtz';

        % The four NV axes in the crystal frame (rows), as in esrFieldSolve.
        NV_AXES = [ 1  1  1; 1 -1 -1; -1  1 -1; -1 -1  1] / sqrt(3);
        NV_THETA = acosd(1/sqrt(3));          % 54.7356 deg, the polar angle of <111>
        GAMMA = 2.8025;                        % MHz/G
        D_DEFAULT = 2870;                      % MHz

        EVENT_OUTPUT_STATE_CHANGED = 'outputStateChanged';
        EVENT_CONTROL_STATE_CHANGED = 'controlStateChanged';
        EVENT_STEP_SIZE_CHANGED = 'stepSizeChanged';
        EVENT_MAGNETIC_FIELD_CHANGED = 'magneticFieldChanged';
        EVENT_CALIBRATION_CHANGED = 'calibrationChanged';
        EVENT_FLIP_STATE_CHANGED = 'flipStateChanged';
        EVENT_PREVIEW_CHANGED = 'previewChanged';

        ESR_WIDTH = 4;          % MHz, half width of the lines in the expected-ESR plot
        I_TOL = 1e-6;           % A: a current this far past a limit counts as at it
    end

    properties (Dependent = true)
        B               % [Bx By Bz], Gauss
        Brho            % |B|, Gauss
        Btheta          % polar angle from +z, degrees
        Bphi            % azimuth from +x, degrees
        Bsphere         % [Brho Btheta Bphi]
        isCalibrated    % true when coilMatrix came from an NV calibration
        frameName       % 'NV crystal frame' or 'nominal coil frame'
        D               % MHz, for predicting lines: the calibration's, else 2870
    end

    properties (SetAccess = private)
        I = [0 0 0];        % [Ix Iy Iz], A -- THE state
        coilMatrix          % M, 3x3, G/A
        coilOffset          % b0, 1x3, G
        calibration = []    % the calibration struct (see setCoilCalibration), [] if nominal
        nominalMatrix       % diag(magnetic_field_ratio) from setupInfo.json
        nominalOffset       % magnetic_field_offset from setupInfo.json
        stateFile = ''      % the json last loaded from or saved to
        previewB = []       % a field the GUI is CONSIDERING (not applied): drawn on the
                            % expected-ESR plot next to the present one. [] = none.
    end

    properties (Access = private)
        psNames
        psChannels
        digitalChannel
        digitalName
        flipAvailable
        flipControl
    end

    properties
        flip            % boolean - flip direction of B (relay; not on setup 7)
        stepSize        % arrow step for B, Gauss and degrees
        stepSizeI       % arrow step for I, A
        maxCurrent      % 1x3, A
        minCurrent      % 1x3, A
        helmOn
        helmControl     % 'ON', 'AUTO' (on during experiments) or 'OFF'
    end

    methods
        function obj = Helmholtz(opts)
            % h = Helmholtz                 the usual: setupInfo.json, supplies, GUI
            % h = Helmholtz(struct('setup', S, 'gui', false))
            %     S replaces the "Helmholtz" section of setupInfo.json (with "type": "dummy"
            %     supplies this runs with no hardware -- see testHelmholtzOffline), and
            %     gui = false skips the window.
            if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
            obj@BaseObject(Helmholtz.NAME);
            obj@EventSender(Helmholtz.NAME);
            obj.helmOn = 0;
            obj.helmControl = obj.CONTROL_STATE{3};
            obj.stepSize = 1;
            obj.stepSizeI = 0.1;
            obj.flip = 0;
            setup = [];
            if isfield(opts, 'setup'); setup = opts.setup; end
            obj.Initialize(setup);
            addBaseObject(obj);  % so it can be reached by getObjByName()
            if isfield(opts, 'gui') && ~opts.gui; return; end
            try
                HelmholtzGUI().start;
            catch err
                obj.close;
                rethrow(err);
            end
        end

        function Initialize(obj, helmholtzStruct)
            if ~exist('helmholtzStruct', 'var') || isempty(helmholtzStruct)
                helmholtzStruct = JsonInfoReader.getJson().Helmholtz;
            end
            obj.psNames = {};
            obj.psChannels = {};
            obj.maxCurrent = zeros(1,3);
            obj.minCurrent = zeros(1,3);
            ratio = zeros(1,3);
            offset = zeros(1,3);
            Inow = zeros(1, 3);
            for k = 1:3
                axisStruct = helmholtzStruct.(obj.AXIS(k));
                psName = sprintf('power_supply_%s', axisStruct.address);
                if ~any(strcmp(obj.psNames, psName))  % one supply may serve several channels
                    powerSupply.create(psName, axisStruct);
                end
                obj.psNames{k} = psName;
                ps = getObjByName(psName);
                if isfield(axisStruct, 'channel')
                    channel = axisStruct.channel;
                else
                    channel = [];
                end
                obj.psChannels{k} = channel;
                if isfield(axisStruct, 'range')
                    ps.setRange(axisStruct.range, channel);
                end
                ps.setMode('current', channel);
                if isfield(axisStruct, 'min_current')
                    minCur = axisStruct.min_current;
                else
                    minCur = 0;
                end
                obj.maxCurrent(k) = min([axisStruct.max_current, ps.maxCurrentActive(channel)]);
                obj.minCurrent(k) = max([minCur, ps.minCurrentActive(channel)]);
                ratio(k) = axisStruct.magnetic_field_ratio;
                offset(k) = axisStruct.magnetic_field_offset;
                Inow(k) = ps.current(ps.chan2ind(channel));
            end
            obj.nominalMatrix = diag(ratio);
            obj.nominalOffset = offset;
            obj.coilMatrix = obj.nominalMatrix;
            obj.coilOffset = obj.nominalOffset;
            obj.calibration = [];

            flipSrtucrt = helmholtzStruct.flip_direction;
            if flipSrtucrt.available
                obj.flipAvailable = 1;
                obj.flipControl = struct;
                obj.flipControl.type = flipSrtucrt.type;
                switch obj.flipControl.type
                    case "power_supply"
                        obj.sendError('power supply control is not implemented');
                    case "pulse_generator"
                        partName = 'flip relay switch';
                        obj.flipControl.switchPhysicalPart = SwitchPgControlled.create(partName, flipSrtucrt.switch);
                        obj.flipControl.switchPhysicalPart.isEnabled = false;
                    case "daq"
                        obj.sendError('DAQ control is not implemented');
                end
            else
                obj.flipAvailable = 0;
            end
            obj.I = Inow;
            obj.setI(obj.maxCurrent*0.5);
        end
    end

    %% the field, from the currents
    methods
        function val = get.B(obj)
            val = (obj.flipSign * obj.coilMatrix * obj.I(:) + obj.coilOffset(:))';
        end

        function val = get.Brho(obj)
            val = norm(obj.B);
        end

        function val = get.Btheta(obj)
            [~, val, ~] = Helmholtz.cart2sphPhys(obj.B);
        end

        function val = get.Bphi(obj)
            [~, ~, val] = Helmholtz.cart2sphPhys(obj.B);
        end

        function val = get.Bsphere(obj)
            [r, t, p] = Helmholtz.cart2sphPhys(obj.B);
            val = [r, t, p];
        end

        function val = get.isCalibrated(obj)
            val = ~isempty(obj.calibration);
        end

        function val = get.frameName(obj)
            if obj.isCalibrated
                val = 'NV crystal frame';
            else
                val = 'nominal coil frame';
            end
        end

        function val = get.D(obj)
            val = obj.D_DEFAULT;
            if ~isempty(obj.calibration) && isfield(obj.calibration, 'D') && ...
                    isfinite(obj.calibration.D)
                val = obj.calibration.D;
            end
        end

        function Iq = currentFor(obj, Bq)
            % The currents that make field Bq (1x3 or Nx3, Gauss). Not range-checked.
            Bq = reshape(Bq, [], 3);
            Iq = (obj.flipSign * (obj.coilMatrix \ (Bq' - obj.coilOffset(:))))';
        end

        function Bq = fieldFor(obj, Iq)
            % The field currents Iq (1x3 or Nx3, A) make.
            Iq = reshape(Iq, [], 3);
            Bq = (obj.flipSign * obj.coilMatrix * Iq' + obj.coilOffset(:))';
        end

        function ok = isReachable(obj, Bq)
            % True for each row of Bq that the supplies can make.
            ok = obj.currentsInRange(obj.currentFor(Bq));
        end

        function set.stepSize(obj, newVal)
            obj.stepSize = newVal;
            obj.sendEvent(struct(obj.EVENT_STEP_SIZE_CHANGED, true));
        end

        function set.stepSizeI(obj, newVal)
            obj.stepSizeI = newVal;
            obj.sendEvent(struct(obj.EVENT_STEP_SIZE_CHANGED, true));
        end
    end

    %% setting the field
    methods
        function setI(obj, I_new)
            I_new = reshape(I_new, 1, 3);
            obj.checkILimits(I_new);
            % within I_TOL of a limit counts as at it: clip, so the supply is never asked
            % for 5.2000001 A because of a rounding error in a field calculation
            I_new = min(max(I_new, obj.minCurrent), obj.maxCurrent);
            obj.I = I_new;
            obj.setCurrent(I_new);
            obj.sendEvent(struct(obj.EVENT_MAGNETIC_FIELD_CHANGED, true));
        end

        function setB(obj, B_new)
            B_new = reshape(B_new, 1, 3);
            Inew = obj.currentFor(B_new);
            if ~obj.currentsInRange(Inew)
                obj.sendError(sprintf(['Out of reach: B = [%.3f %.3f %.3f] G (%s) needs ' ...
                    'I = [%.3f %.3f %.3f] A, and the supplies give [%.2f %.2f %.2f] to ' ...
                    '[%.2f %.2f %.2f] A.'], B_new, obj.frameName, Inew, obj.minCurrent, ...
                    obj.maxCurrent));
            end
            obj.setI(Inew);
        end

        function setBsphere(obj, rho, theta, phi)
            % Brho in Gauss, Btheta (polar, from z) and Bphi (azimuth) in degrees.
            obj.checkTheta(theta);
            obj.setB(Helmholtz.sph2cartPhys(rho, theta, phi));
        end

        function setBrho(obj, Brho_new)
            if Brho_new < 0
                obj.sendError('Brho must be >= 0. To reverse the field, add 180 to Bphi and take 180 - Btheta.');
            end
            [~, t, p] = obj.anglesForSetting();
            obj.setB(Helmholtz.sph2cartPhys(Brho_new, t, p));
        end

        function setBtheta(obj, Btheta_new)
            obj.checkTheta(Btheta_new);
            [r, ~, p] = obj.anglesForSetting();
            obj.setB(Helmholtz.sph2cartPhys(r, Btheta_new, p));
        end

        function setBphi(obj, Bphi_new)
            [r, t, ~] = obj.anglesForSetting();
            obj.setB(Helmholtz.sph2cartPhys(r, t, Bphi_new));
        end

        function relativeBsphereStep(obj, index, step)
            switch index
                case 1
                    obj.set('Brho', obj.Brho + step);
                case 2
                    obj.set('Btheta', obj.Btheta + step);
                case 3
                    obj.set('Bphi', obj.Bphi + step);
            end
        end

        function relativeBcartStep(obj, index, step)
            Bn = obj.B;
            Bn(index) = Bn(index) + step;
            obj.setB(Bn);
        end

        function relativeCurrentStep(obj, index, step)
            In = obj.I;
            In(index) = In(index) + step;
            obj.setI(In);
        end

        function set(obj, what, value)
            what = char(what);   % "Bx" (string) as well as 'Bx'
            switch what
                case {'Bx', 'By', 'Bz'}
                    Bn = obj.B;
                    Bn(what(2) - 'w') = value;          % 'x'->1, 'y'->2, 'z'->3
                    obj.setB(Bn);
                case {'Ix', 'Iy', 'Iz'}
                    In = obj.I;
                    In(what(2) - 'w') = value;
                    obj.setI(In);
                case 'Brho'
                    obj.setBrho(value);
                case 'Btheta'
                    obj.setBtheta(value);
                case 'Bphi'
                    obj.setBphi(value);
                otherwise
                    obj.sendError(sprintf('no field %s to set', what));
            end
        end

        function Bq = settingToField(obj, what, value)
            what = char(what);   % "Bx" (string) as well as 'Bx'
            % The field that set(what, value) would ask for, without setting it.
            [r, t, p] = obj.anglesForSetting();
            switch what
                case {'Bx', 'By', 'Bz'}
                    Bq = obj.B; Bq(what(2) - 'w') = value;
                case {'Ix', 'Iy', 'Iz'}
                    In = obj.I; In(what(2) - 'w') = value; Bq = obj.fieldFor(In);
                case 'Brho'
                    Bq = Helmholtz.sph2cartPhys(value, t, p);
                case 'Btheta'
                    Bq = Helmholtz.sph2cartPhys(r, value, p);
                case 'Bphi'
                    Bq = Helmholtz.sph2cartPhys(r, t, value);
                otherwise
                    error('Helmholtz:what', 'no field %s', what);
            end
        end
    end

    %% the NV axes
    methods
        function [lo, hi] = nvReach(obj, k)
            % How far along signed NV direction k (1..8, see nvDirections) the coils reach:
            % the field can be r*dir for r in [lo, hi] Gauss, or [NaN NaN] if never.
            dirs = Helmholtz.nvDirections();
            u = dirs(k, :);
            [lo, hi] = obj.rayLimits(u);
        end

        function alignToNV(obj, k, rho)
            % Put B along an NV axis. k: 1..8 (Helmholtz.nvDirections) or a direction such as
            % [1 1 1] or [-1 1 1]. rho: |B| in Gauss (default: the present |B|).
            %   h.alignToNV([1 1 1], 30)     % 30 G along [111]: theta 54.7356, phi 45
            if ~exist('rho', 'var') || isempty(rho); rho = obj.Brho; end
            u = obj.nvDirectionFrom(k);
            if ~obj.isCalibrated
                obj.sendWarning(['No coil calibration loaded: B is in the NOMINAL coil frame, ' ...
                    'so this is only roughly along an NV. Run calibCoilMatrix (level 3) first.']);
            end
            [lo, hi] = obj.rayLimits(u);
            if ~(rho >= lo - 1e-9 && rho <= hi + 1e-9)
                obj.sendError(sprintf(['|B| = %.2f G along [%s] is out of reach: the supplies ' ...
                    'allow %.2f to %.2f G in that direction.'], rho, ...
                    num2str(round(u*sqrt(3))), lo, hi));
            end
            obj.setB(rho * u);
        end

        function [k, cosAng] = nearestNV(obj)
            % The signed NV direction (1..8) closest to the present field, and cos of the angle.
            dirs = Helmholtz.nvDirections();
            c = dirs * obj.B(:) / max(obj.Brho, eps);
            [cosAng, k] = max(c);
        end

        function BtoNV(obj)
            % Turn the field onto the nearest NV axis at the present |B|. If that |B| is not
            % reachable along it, the nearest reachable |B| is used, and said.
            k = obj.nearestNV();
            [lo, hi] = obj.nvReach(k);
            rho = obj.Brho;
            if isnan(lo)
                obj.sendError('The nearest NV axis cannot be reached at any |B| with these supplies.');
            end
            if rho < lo || rho > hi
                rhoNew = min(max(rho, lo), hi);
                obj.sendWarning(sprintf(['|B| = %.2f G is out of reach along that axis ' ...
                    '(%.2f to %.2f G); using %.2f G.'], rho, lo, hi, rhoNew));
                rho = rhoNew;
            end
            obj.alignToNV(k, rho);
        end

        function T = nvTable(obj)
            % One row per signed NV direction: its index, direction, theta, phi, the |B|
            % range the coils reach along it, and the angle to the present field.
            dirs = Helmholtz.nvDirections();
            T = struct('k', {}, 'dir', {}, 'theta', {}, 'phi', {}, 'rhoMin', {}, ...
                       'rhoMax', {}, 'angleToB', {});
            for k = 1:8
                [~, t, p] = Helmholtz.cart2sphPhys(dirs(k, :));
                [lo, hi] = obj.nvReach(k);
                ang = acosd(max(min(dirs(k, :) * obj.B(:) / max(obj.Brho, eps), 1), -1));
                T(k) = struct('k', k, 'dir', round(dirs(k, :)*sqrt(3)), 'theta', t, ...
                              'phi', p, 'rhoMin', lo, 'rhoMax', hi, 'angleToB', ang);
            end
        end

        function f = expectedLines(obj, Bq)
            % The 8 NV transition frequencies (MHz) at field Bq (default: the present one).
            % Exact spin-1 Hamiltonian, the same physics as esrFieldSolve. Only means
            % something in the NV crystal frame, i.e. after a coil calibration.
            if ~exist('Bq', 'var'); Bq = obj.B; end
            f = Helmholtz.nvLines(Bq, obj.D);
        end

        % ---- targets: NV axes, and the cube axes where all four NVs are equal ----------

        function T = targetTable(obj)
            % Every direction the GUI offers (Helmholtz.targetDirections), with the |B| range
            % the supplies reach EXACTLY along it ([NaN NaN] if never) and the angle to the
            % present field.
            T = Helmholtz.targetDirections();
            for k = 1:numel(T)
                [lo, hi] = obj.rayLimits(T(k).u);
                T(k).rhoMin = lo;
                T(k).rhoMax = hi;
                T(k).angleToB = acosd(max(min(T(k).u * obj.B(:) / max(obj.Brho, eps), 1), -1));
            end
        end

        function [Bq, Iq, angDeg, exact] = closestReachable(obj, Btarget)
            % The reachable field closest to Btarget (Euclidean, in Gauss), the currents for
            % it, its angle to Btarget, and whether it IS Btarget.
            %
            % Along a cube axis the supplies often need a few mA of NEGATIVE current to cancel
            % the zero-current field and the cross-talk exactly (setup 7: -0.04 A of Ix for
            % 40 G along z). Unipolar supplies cannot, but the closest reachable field is then
            % a fraction of a degree off -- far less than a linewidth of splitting -- and that
            % is what this returns. A box-constrained least squares in three currents, solved
            % exactly by trying every set of active limits (27) and keeping the best feasible.
            Btarget = reshape(Btarget, 3, 1);
            M = obj.flipSign * obj.coilMatrix;
            r = Btarget - obj.coilOffset(:);
            lo = obj.minCurrent(:); hi = obj.maxCurrent(:);
            best = inf; Iq = [];
            for code = 0:26
                s = mod(floor(code ./ [1 3 9]), 3);      % 0 free, 1 at min, 2 at max
                Ic = zeros(3, 1);
                fixed = s > 0;
                Ic(s == 1) = lo(s == 1); Ic(s == 2) = hi(s == 2);
                if any(~fixed)
                    A = M(:, ~fixed);
                    Ic(~fixed) = A \ (r - M(:, fixed) * Ic(fixed));
                end
                if any(Ic < lo - obj.I_TOL | Ic > hi + obj.I_TOL); continue; end
                d = norm(M * Ic - r);
                if d < best - 1e-12; best = d; Iq = Ic; end
            end
            Iq = min(max(Iq', obj.minCurrent), obj.maxCurrent);
            Bq = obj.fieldFor(Iq);
            angDeg = acosd(max(min(Bq * Btarget / max(norm(Bq) * norm(Btarget), eps), 1), -1));
            exact = best < 1e-6;
        end

        function setDirection(obj, u, rho, maxDeg)
            % B = rho*u if the supplies can make it, else the closest field they can, as long
            % as it is within maxDeg degrees of u (default 1). Prints the miss.
            %   h.setDirection([0 0 1], 30)        % along +z: all four NVs equal
            if ~exist('maxDeg', 'var') || isempty(maxDeg); maxDeg = 1; end
            u = reshape(u, 1, 3) / norm(u);
            [Bq, Iq, ang, exact] = obj.closestReachable(rho * u);
            if ~exact
                if ang > maxDeg
                    obj.sendError(sprintf(['%.2f G along [%s] is out of reach, and the closest ' ...
                        'reachable field is %.2f deg off (limit %.2f).'], rho, ...
                        num2str(u, '%.3f '), ang, maxDeg));
                end
                fprintf(['Helmholtz: %.2f G along [%s] needs currents outside the supplies'' ' ...
                         'range; set the closest reachable field instead: %.3f G, %.3f deg off.\n'], ...
                        rho, num2str(u, '%.3f '), norm(Bq), ang);
            end
            obj.setI(Iq);
        end

        function setPreview(obj, Bq)
            % A field being considered, drawn next to the present one. [] clears it.
            if ~isempty(Bq); Bq = reshape(Bq, 1, 3); end
            obj.previewB = Bq;
            obj.sendEvent(struct(obj.EVENT_PREVIEW_CHANGED, true));
        end

        function flipDirection(obj, state)
            if obj.flipAvailable
                state = num2str(state);
                if ~any(strcmpi({'on', 'off', '0', '1'}, state))
                    obj.sendError('flip state can get just ''ON'' or ''OFF''')
                end
                obj.flip = BooleanHelper.onOffToBool(state);
                switch obj.flipControl.type
                    case "power_supply"
                        obj.sendError('power supply control is not implemented');
                    case "pulse_generator"
                        obj.flipControl.switchPhysicalPart.isEnabled = obj.flip;
                    case "daq"
                        obj.sendError('DAQ control is not implemented');
                end
                obj.sendEvent(struct(obj.EVENT_FLIP_STATE_CHANGED, true));
                obj.sendEvent(struct(obj.EVENT_MAGNETIC_FIELD_CHANGED, true));
            else
                obj.sendWarning('Flip direction option is not available');
            end
        end
    end

    %% the coil calibration
    methods
        function setCoilCalibration(obj, cal)
            % Use a coil calibration: cal is coilMatrixFit's out.calibration (or the
            % "calibration" block of the json), with at least
            %   matrix  3x3 G/A     offset  1x3 G
            % and optionally date, D, rmsMHz, matrixErr, offsetErr, nSpectra, currentMin/Max.
            % The CURRENTS are kept, so the coils do not move -- what changes is what we know
            % the field to be, and the frame it is expressed in.
            if isfield(cal, 'calibration') && ~isfield(cal, 'matrix')
                cal = cal.calibration;          % accept coilMatrixFit's whole output
            end
            Mc = cal.matrix;
            if ~isequal(size(Mc), [3 3]) || ~all(isfinite(Mc(:))) || rcond(Mc) < 1e-6
                obj.sendError('coil calibration: matrix must be a finite, invertible 3x3 (G/A).');
            end
            b = reshape(cal.offset, 1, 3);
            if ~all(isfinite(b))
                obj.sendError('coil calibration: offset must be 3 finite numbers (G).');
            end
            cal.matrix = Mc;
            cal.offset = b;
            obj.coilMatrix = Mc;
            obj.coilOffset = b;
            obj.calibration = cal;
            obj.sendEvent(struct(obj.EVENT_CALIBRATION_CHANGED, true));
            obj.sendEvent(struct(obj.EVENT_MAGNETIC_FIELD_CHANGED, true));
        end

        function clearCoilCalibration(obj)
            % Back to the nominal model from setupInfo.json. The currents are kept.
            obj.coilMatrix = obj.nominalMatrix;
            obj.coilOffset = obj.nominalOffset;
            obj.calibration = [];
            obj.sendEvent(struct(obj.EVENT_CALIBRATION_CHANGED, true));
            obj.sendEvent(struct(obj.EVENT_MAGNETIC_FIELD_CHANGED, true));
        end

        function s = calibrationSummary(obj)
            % One line for the GUI and the log.
            if ~obj.isCalibrated
                s = 'UNCALIBRATED: nominal coil model from setupInfo.json (coil frame)';
                return
            end
            c = obj.calibration;
            s = 'calibrated';
            if isfield(c, 'date') && ~isempty(c.date); s = [s ' ' char(c.date)]; end
            if isfield(c, 'rmsMHz') && isfinite(c.rmsMHz)
                s = sprintf('%s, line rms %.2f MHz', s, c.rmsMHz);
            end
            if isfield(c, 'nSpectra'); s = sprintf('%s, %d spectra', s, c.nSpectra); end
            s = [s ' (NV crystal frame)'];
        end

        function w = extrapolationWarning(obj)
            % Non-empty when the present currents are outside what the calibration measured.
            w = '';
            if ~obj.isCalibrated || ~isfield(obj.calibration, 'currentMin'); return; end
            lo = reshape(obj.calibration.currentMin, 1, 3);
            hi = reshape(obj.calibration.currentMax, 1, 3);
            out = obj.I < lo - 0.05 | obj.I > hi + 0.05;
            if any(out)
                w = sprintf(['coil %s current outside the calibrated range -- the model is ' ...
                    'linear, so this is usually fine, but it is an extrapolation'], ...
                    obj.AXIS(out));
            end
        end
    end

    %% limits
    methods
        function [lo, hi] = limits(obj, what)
            what = char(what);   % "Bx" (string) as well as 'Bx'
            % The range of one coordinate with the other two held at their present values.
            % what: 'Bx' 'By' 'Bz' 'Brho' 'Btheta' 'Bphi' 'Ix' 'Iy' 'Iz'. [NaN NaN] if the
            % present state is itself out of range. For Bphi the interval may run past
            % +-180 (it is the stretch of the circle that contains the present phi).
            switch what
                case {'Ix', 'Iy', 'Iz'}
                    k = what(2) - 'w';
                    lo = obj.minCurrent(k); hi = obj.maxCurrent(k);
                otherwise
                    [lo, hi] = Helmholtz.fieldLimits(obj.flipSign * obj.coilMatrix, ...
                        obj.coilOffset, obj.minCurrent, obj.maxCurrent, obj.B, what);
            end
        end

        function lim = BsphereLowerLimits(obj)
            lim = [obj.limitLo('Brho'), obj.limitLo('Btheta'), obj.limitLo('Bphi')];
        end

        function lim = BsphereUpperLimits(obj)
            lim = [obj.limitHi('Brho'), obj.limitHi('Btheta'), obj.limitHi('Bphi')];
        end

        function lim = BcartLowerLimits(obj)
            lim = [obj.limitLo('Bx'), obj.limitLo('By'), obj.limitLo('Bz')];
        end

        function lim = BcartUpperLimits(obj)
            lim = [obj.limitHi('Bx'), obj.limitHi('By'), obj.limitHi('Bz')];
        end

        function OK = checkILimits(obj, I, dispError)
            if ~exist('dispError', 'var'); dispError = 1; end
            OK = obj.currentsInRange(I);
            if ~OK
                err = sprintf(['Out of limits! The current range is [%.2f %.2f %.2f] to ' ...
                    '[%.2f %.2f %.2f] A; asked for [%.3f %.3f %.3f] A.'], obj.minCurrent, ...
                    obj.maxCurrent, reshape(I(1, :), 1, 3));
                obj.dispErrorWarning(err, dispError);
            end
        end

        function OK = checkLimits(obj, what, value, dispError)
            what = char(what);   % "Bx" (string) as well as 'Bx'
            % True when EVERY value in `value` can be set for coordinate `what` (the other
            % two held). For 'I' and 'B', value is a 1x3 or Nx3 list of whole settings.
            if ~exist('dispError', 'var'); dispError = 1; end
            switch what
                case 'I'
                    OK = all(obj.currentsInRange(reshape(value, [], 3)));
                case 'B'
                    OK = all(obj.isReachable(reshape(value, [], 3)));
                case {'Bx', 'By', 'Bz', 'Brho', 'Btheta', 'Bphi', 'Ix', 'Iy', 'Iz'}
                    OK = true;
                    for v = value(:)'
                        if strcmp(what, 'Btheta') && (v < 0 || v > 180); OK = false; break; end
                        if strcmp(what, 'Brho') && v < 0; OK = false; break; end
                        if ~obj.isReachable(obj.settingToField(what, v)); OK = false; break; end
                    end
                otherwise
                    obj.sendError(sprintf('no type of field %s', what));
            end
            if ~OK
                if any(strcmp(what, {'I', 'B'}))
                    err = sprintf('Out of limits! Some %s settings cannot be reached.', what);
                else
                    [lo, hi] = obj.limits(what);
                    err = sprintf(['Out of limits! With the other coordinates held, %s can go ' ...
                        'from %.3f to %.3f.'], what, lo, hi);
                end
                obj.dispErrorWarning(err, dispError);
            end
        end

        function dispErrorWarning(obj, err, dispError)
            % 'err': error string;
            % 'dispError': 0 - no display, 1 - error, 2 - warning.
            if dispError == 1
                obj.sendError(err)
            elseif dispError == 2
                warning('Helmholtz:limits', '%s', err)
            end
        end
    end

    %% communication to power supplies methods
    methods
        function setCurrent(obj, current)
            for k = 1:3
                ps = getObjByName(obj.psNames{k});
                channel = obj.psChannels{k};
                pause(0.1);
                ps.setCurrent(current(k), channel);
            end
        end

        function output(obj, value)
            for k = 1:3
                ps = getObjByName(obj.psNames{k});
                channel = obj.psChannels{k};
                ps.output(value, channel);
            end
            obj.helmOn = strcmpi(value, 'on') || any((value==1));
            obj.sendEvent(struct(obj.EVENT_OUTPUT_STATE_CHANGED, true));
        end

        function control(obj, value)
            % 'ON', 'AUTO' (on only while an experiment runs) or 'OFF'. Case-insensitive --
            % it used to validate case-insensitively and then branch case-SENSITIVELY, so
            % control('on') passed and switched the supplies OFF.
            value = upper(char(value));
            if ~any(strcmp(value, obj.CONTROL_STATE))
                error('control state can be just ''ON'', ''AUTO'' or ''OFF''')
            end
            obj.helmControl = value;
            if strcmp(value, obj.CONTROL_STATE{1})
                obj.output('on');
            else
                obj.output('off');
            end
            obj.sendEvent(struct(obj.EVENT_CONTROL_STATE_CHANGED, true));
        end

        function close(obj)
            for k = 1:numel(obj.psNames)
                ps = getObjByName(obj.psNames{k});
                if ~isempty(ps)
                    ps.close;
                end
            end
            h = getObjByName('Helmholtz');
            if ~isempty(h); h.delete; end
        end

        function reset(obj)
            for k = 1:3
                ps = getObjByName(obj.psNames{k});
                ps.reset;
            end
            obj.setI(max(obj.minCurrent, [0 0 0]));
        end
    end

    %% save / load -- the "Helmholtz" block of expProperties.json
    methods
        function s = stateStruct(obj)
            % What goes into the json. Plain numbers, strings and a nested calibration.
            s = struct();
            s.I = obj.I;
            s.B = obj.B;
            s.frame = obj.frameName;
            s.Brho = obj.Brho;
            s.Btheta = obj.Btheta;
            s.Bphi = obj.Bphi;
            s.flip = logical(obj.flip);
            s.stepSize = obj.stepSize;
            s.stepSizeI = obj.stepSizeI;
            if obj.isCalibrated
                c = Helmholtz.cleanCalibration(obj.calibration);
                % json has no NaN: drop a non-finite field (an error estimate, say) rather
                % than refuse to save the calibration itself
                for f = fieldnames(c)'
                    v = c.(f{1});
                    if isnumeric(v) && any(~isfinite(v(:)))
                        if any(strcmp(f{1}, {'matrix', 'offset'}))
                            obj.sendError('Helmholtz.save: the coil calibration itself is not finite.');
                        end
                        fprintf('Helmholtz.save: calibration.%s is not finite; not saved.\n', f{1});
                        c = rmfield(c, f{1});
                    end
                end
                s.calibration = c;
            end
        end

        function save(obj, address)
            % h.save(jsonFile): write the state and the coil calibration into the
            % "Helmholtz" block of that expProperties.json (created if absent), through
            % expPropsWrite: backed up first, verified after, nothing else touched.
            if ~exist('address', 'var') || isempty(address); address = obj.stateFile; end
            if isempty(address)
                obj.sendError('Helmholtz.save: no file given, and none was loaded earlier.');
            end
            [~, ~, ext] = fileparts(address);
            if ~strcmpi(ext, '.json')
                obj.sendError(['Helmholtz.save writes to expProperties.json now (key "Helmholtz"); ' ...
                    'helmholtzState.mat is no longer written. Pass the json, e.g. setup7Paths().json.']);
            end
            if exist('expPropsWrite', 'file') ~= 2
                obj.sendError(['Helmholtz.save needs expPropsWrite (Setup 7\functions) on the ' ...
                    'path -- run setup7AddPath.']);
            end
            if ~isfile(address)
                fid = fopen(address, 'w');
                if fid < 0; obj.sendError(sprintf('cannot create %s', address)); end
                fprintf(fid, '{\n}\n');
                fclose(fid);
            end
            res = expPropsWrite(struct(obj.JSON_KEY, obj.stateStruct()), ...
                struct('jsonLocation', address, 'create', true, 'verbose', false));
            if ~res.ok
                obj.sendError(sprintf('Helmholtz.save: the json was not written (%s).', address));
            end
            obj.stateFile = address;
            fprintf('Helmholtz state saved to %s ("%s"): I = [%.4f %.4f %.4f] A, B = [%.3f %.3f %.3f] G (%s)\n', ...
                    address, obj.JSON_KEY, obj.I, obj.B, obj.frameName);
        end

        function load(obj, address)
            % h.load(jsonFile): the "Helmholtz" block of an expProperties.json -- the coil
            % calibration if it has one, then the currents.
            % h.load(matFile): a legacy helmholtzState.mat. Its B was in the NOMINAL coil
            % frame, so it is turned back into the currents that made it; any calibration
            % already loaded is kept.
            [~, ~, ext] = fileparts(address);
            if strcmpi(ext, '.mat')
                obj.loadLegacyMat(address);
                return
            end
            J = jsondecode(fileread(address));
            if ~isfield(J, obj.JSON_KEY)
                % An older results folder: the state is still in helmholtzState.mat next to
                % the json. Load that, and say how to move it into the json.
                legacy = fullfile(fileparts(address), 'helmholtzState.mat');
                if isfile(legacy)
                    fprintf('Helmholtz.load: %s has no "%s" block; using %s.\n', address, ...
                            obj.JSON_KEY, legacy);
                    obj.loadLegacyMat(legacy);
                    obj.stateFile = address;    % so the next Save moves it into the json
                    return
                end
                obj.sendError(sprintf(['%s has no "%s" block (and no helmholtzState.mat beside ' ...
                    'it). Set the field and save it with h.save(''%s'').'], address, ...
                    obj.JSON_KEY, address));
            end
            s = J.(obj.JSON_KEY);
            % first, so the window's redraws below already see the new file (the
            % expected-ESR plot reads the operating frequency from it)
            obj.stateFile = address;
            obj.previewB = [];
            if isfield(s, 'calibration') && isstruct(s.calibration)
                obj.setCoilCalibration(Helmholtz.cleanCalibration(s.calibration));
            else
                obj.clearCoilCalibration();
            end
            if isfield(s, 'flip') && obj.flipAvailable
                obj.flipDirection(double(s.flip));
            end
            if isfield(s, 'stepSize'); obj.stepSize = s.stepSize; end
            if isfield(s, 'stepSizeI'); obj.stepSizeI = s.stepSizeI; end
            if isfield(s, 'I')
                obj.setI(reshape(s.I, 1, 3));
                if isfield(s, 'B')
                    dB = norm(obj.B - reshape(s.B, 1, 3));
                    if dB > 0.01
                        fprintf(['Helmholtz.load: the saved B differs from what the saved currents ' ...
                                 'make with this calibration by %.3f G -- the calibration changed ' ...
                                 'after the state was saved. The CURRENTS were restored.\n'], dB);
                    end
                end
            elseif isfield(s, 'B')
                obj.setB(reshape(s.B, 1, 3));
            end
            obj.stateFile = address;
            fprintf('Helmholtz state loaded from %s: I = [%.4f %.4f %.4f] A, |B| = %.3f G, theta %.3f, phi %.3f deg (%s)\n', ...
                    address, obj.I, obj.Brho, obj.Btheta, obj.Bphi, obj.frameName);
        end

        function loadLegacyMat(obj, address)
            % helmholtzState.mat: expStruct.B in the nominal coil frame, plus flip and steps.
            S = load(address);
            e = S.expStruct;
            Iold = (reshape(e.B, 1, 3) - obj.nominalOffset) ./ diag(obj.nominalMatrix)';
            if isfield(e, 'flip') && obj.flipAvailable; obj.flipDirection(e.flip); end
            if isfield(e, 'stepSize'); obj.stepSize = e.stepSize; end
            if isfield(e, 'stepSizeI'); obj.stepSizeI = e.stepSizeI; end
            obj.setI(Iold);
            fprintf(['Helmholtz: legacy %s loaded -- its B = [%.3f %.3f %.3f] G was in the nominal ' ...
                     'coil frame, i.e. I = [%.4f %.4f %.4f] A. To move it into the json: ' ...
                     'h.save(setup7Paths().json)\n'], address, e.B, obj.I);
        end
    end

    %% plots
    methods
        function xl = plotExpectedESR(obj, gAxes, opts)
            % The ESR this field should give, from the exact Hamiltonian -- and, when the GUI
            % is considering a new field (previewB), that one too, dashed. Every dip is marked
            % and labelled with its frequency and how many orientations share it.
            %   opts.operatingMHz   the json "frequency", drawn as the MW line (optional)
            % Returns the x limits it chose, so a caller can keep a user's zoom.
            if ~exist('gAxes', 'var') || isempty(gAxes)
                fig = figure('Color', [252 252 251]/255, 'Name', 'expected ESR');
                gAxes = axes('Parent', fig);
            end
            if ~exist('opts', 'var') || isempty(opts); opts = struct(); end
            C = Helmholtz.plotColours();
            ax = gAxes;
            % NOT cla: cla leaves objects with hidden handles (the fill, the dip markers),
            % so every redraw would stack on the last
            delete(allchild(ax));
            legend(ax, 'off');
            hold(ax, 'on');
            w = obj.ESR_WIDTH;
            D0 = obj.D;
            series = {struct('B', obj.B, 'col', C.s1, 'name', 'present field', 'style', '-')};
            if ~isempty(obj.previewB)
                series{end+1} = struct('B', obj.previewB, 'col', C.s2, ...
                    'name', 'preview (not applied)', 'style', '--');
            end
            allDips = [];
            for s = 1:numel(series)
                L = Helmholtz.nvLinesByAxis(series{s}.B, D0);
                series{s}.dips = Helmholtz.groupDips(L(:));
                allDips = [allDips; series{s}.dips(:, 1)];                    %#ok<AGROW>
            end
            lo = min([allDips; D0]) - 7*w; hi = max([allDips; D0]) + 7*w;
            if isfield(opts, 'operatingMHz') && isfinite(opts.operatingMHz)
                lo = min(lo, opts.operatingMHz - 3*w); hi = max(hi, opts.operatingMHz + 3*w);
            end
            f = linspace(lo, hi, 3000);
            depth1 = 0.01;
            yMin = 1;
            H = gobjects(0); names = {};
            for s = 1:numel(series)
                dp = series{s}.dips;
                y = ones(size(f));
                for k = 1:size(dp, 1)
                    y = y - depth1 * dp(k, 2) * w^2 ./ ((f - dp(k, 1)).^2 + w^2);
                end
                yMin = min(yMin, min(y));
                if s == 1
                    fill(ax, [f, fliplr(f)], [y, ones(size(f))], series{s}.col, ...
                         'FaceAlpha', 0.14, 'EdgeColor', 'none', 'HandleVisibility', 'off', ...
                         'HitTest', 'off');
                end
                hl = plot(ax, f, y, series{s}.style, 'Color', series{s}.col, 'LineWidth', 2);
                H(end+1) = hl; names{end+1} = series{s}.name;                  %#ok<AGROW>
                series{s}.y = y;
            end
            span = 1 - yMin; if span <= 0; span = 0.01; end
            % the dips: a marker at each minimum, and its frequency (x multiplicity) under it
            for s = 1:numel(series)
                dp = series{s}.dips;
                yd = interp1(f, series{s}.y, dp(:, 1));
                hm = plot(ax, dp(:, 1), yd, 'v', 'MarkerSize', 7, 'LineWidth', 1, ...
                          'Color', series{s}.col, 'MarkerFaceColor', series{s}.col, ...
                          'HandleVisibility', 'off');
                try
                    hm.DataTipTemplate.DataTipRows(1).Label = 'dip (MHz)';
                    hm.DataTipTemplate.DataTipRows(2).Label = 'signal';
                    hm.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('orientations', dp(:, 2));
                catch
                end
                for k = 1:size(dp, 1)
                    lbl = sprintf('%.1f', dp(k, 1));
                    if dp(k, 2) > 1; lbl = sprintf('%s  x%d', lbl, dp(k, 2)); end
                    if s == 1
                        % present field: under its dip
                        text(ax, dp(k, 1), yd(k) - 0.04*span, lbl, 'Color', series{s}.col, ...
                             'FontSize', 8, 'FontWeight', 'bold', 'Rotation', 90, ...
                             'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                             'HitTest', 'off');
                    else
                        % preview: above the baseline, so a close pair of labels cannot overlap
                        text(ax, dp(k, 1), 1 + 0.04*span, lbl, 'Color', series{s}.col, ...
                             'FontSize', 8, 'FontWeight', 'bold', 'Rotation', 90, ...
                             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                             'HitTest', 'off');
                    end
                end
            end
            top = 0.12; if numel(series) > 1; top = 0.5; end
            yl = [1 - span*1.75, 1 + span*top];
            plot(ax, D0*[1 1], yl, ':', 'Color', C.ink2, 'LineWidth', 1, 'HandleVisibility', 'off', ...
                 'HitTest', 'off');
            text(ax, D0, yl(1), sprintf(' D %.1f', D0), 'Color', C.ink2, 'FontSize', 8, ...
                 'VerticalAlignment', 'bottom', 'HitTest', 'off');
            if isfield(opts, 'operatingMHz') && isfinite(opts.operatingMHz)
                ho = plot(ax, opts.operatingMHz*[1 1], yl, '-.', 'Color', C.ok, 'LineWidth', 1.3);
                text(ax, opts.operatingMHz, yl(1), sprintf(' MW %.2f', opts.operatingMHz), ...
                     'Color', C.ok, 'FontSize', 8, 'FontWeight', 'bold', ...
                     'VerticalAlignment', 'bottom', 'HitTest', 'off');
                H(end+1) = ho; names{end+1} = 'operating frequency (json)';
            end
            hold(ax, 'off');
            Helmholtz.styleAxes(ax, C);
            xlim(ax, [lo hi]); ylim(ax, yl);
            xl = [lo hi];
            xlabel(ax, 'frequency (MHz)', 'Color', C.ink2);
            ylabel(ax, 'expected signal (normalised)', 'Color', C.ink2);
            [kN, c] = obj.nearestNV();
            dirs = Helmholtz.nvDirections();
            t1 = sprintf('|B| = %.2f G,  %.2f deg from NV [%+d %+d %+d]', obj.Brho, ...
                         acosd(min(c, 1)), round(dirs(kN, :)*sqrt(3)));
            if ~obj.isCalibrated
                t1 = [t1 '   --   UNCALIBRATED: lines assume the coil frame is the crystal frame'];
            end
            title(ax, t1, 'Color', C.ink, 'FontWeight', 'bold', 'FontSize', 10);
            % always set: a subtitle is not a child, so it outlives the redraw
            if ~isempty(obj.previewB)
                subtitle(ax, sprintf('preview: |B| = %.2f G  (press Go to apply)', ...
                         norm(obj.previewB)), 'Color', C.s2, 'FontSize', 9);
            else
                subtitle(ax, '');
            end
            lg = legend(ax, H, names, 'Location', 'southeast', 'Box', 'off', 'FontSize', 8);
            lg.TextColor = C.ink2;
        end
    end

    %% helpers
    methods (Access = private)
        function s = flipSign(obj)
            s = 1;
            if ~isempty(obj.flipAvailable) && obj.flipAvailable && obj.flip; s = -1; end
        end

        function ok = currentsInRange(obj, Iq)
            Iq = reshape(Iq, [], 3);
            ok = all(Iq >= obj.minCurrent - obj.I_TOL & Iq <= obj.maxCurrent + obj.I_TOL, 2);
        end

        function [r, t, p] = anglesForSetting(obj)
            % The present (rho, theta, phi), with a usable direction even at |B| = 0.
            [r, t, p] = Helmholtz.cart2sphPhys(obj.B);
        end

        function checkTheta(obj, theta)
            if ~(theta >= 0 && theta <= 180)
                obj.sendError(sprintf(['Btheta is the POLAR angle from +z and must be within ' ...
                    '0..180 deg (asked for %.3f). The azimuth is Bphi.'], theta));
            end
        end

        function [lo, hi] = rayLimits(obj, u)
            [lo, hi] = Helmholtz.rayInterval(obj.flipSign * obj.coilMatrix, obj.coilOffset, ...
                obj.minCurrent, obj.maxCurrent, u);
        end

        function u = nvDirectionFrom(obj, k)
            dirs = Helmholtz.nvDirections();
            if isscalar(k)
                if k < 1 || k > 8 || k ~= round(k)
                    obj.sendError('NV direction index must be 1..8 (see Helmholtz.nvDirections).');
                end
                u = dirs(k, :);
            else
                v = reshape(k, 1, 3);
                if ~all(abs(abs(v) - 1) < 1e-9)
                    obj.sendError('NV direction must be like [1 1 1] or [-1 1 1] (all entries +-1).');
                end
                u = v / sqrt(3);
            end
        end

        function c = nearestNVcos(obj)
            [~, c] = obj.nearestNV();
        end

        function v = limitLo(obj, what)
            [v, ~] = obj.limits(what);
        end

        function v = limitHi(obj, what)
            [~, v] = obj.limits(what);
        end
    end

    methods (Static)
        function dirs = nvDirections()
            % The eight signed NV directions, crystal frame, unit rows. 1-4 point up (+z,
            % theta = 54.7356, phi = 45, 135, -135, -45); 5-8 are their opposites.
            v = [1 1 1; -1 1 1; -1 -1 1; 1 -1 1];
            dirs = [v; -v] / sqrt(3);
        end

        function T = targetDirections()
            % The directions the GUI offers to put B along:
            %   'NV'    the 8 signed NV axes -- one orientation aligned, three equivalent at
            %           109.47 deg: four lines, the outer pair split by 2*gamma*|B|
            %   'equal' the 6 cube axes (+-x, +-y, +-z) -- every NV axis at 54.74 deg to B,
            %           so all four orientations are equivalent: TWO lines, 4x the depth
            dirs = Helmholtz.nvDirections();
            T = struct('label', {}, 'kind', {}, 'u', {}, 'rhoMin', {}, 'rhoMax', {}, ...
                       'angleToB', {});
            for k = 1:8
                T(end+1) = struct('label', sprintf('NV [%+d %+d %+d]', round(dirs(k, :)*sqrt(3))), ...
                    'kind', 'NV', 'u', dirs(k, :), 'rhoMin', NaN, 'rhoMax', NaN, ...
                    'angleToB', NaN);                                                    %#ok<AGROW>
            end
            names = {'+x', '+y', '+z', '-x', '-y', '-z'};
            E = [eye(3); -eye(3)];
            for k = 1:6
                T(end+1) = struct('label', sprintf('%s  (all 4 NV equal)', names{k}), ...
                    'kind', 'equal', 'u', E(k, :), 'rhoMin', NaN, 'rhoMax', NaN, ...
                    'angleToB', NaN);                                                    %#ok<AGROW>
            end
        end

        function L = nvLinesByAxis(Bq, D)
            % 4 x 2: [lower upper] transition (MHz) of each NV orientation, in NV_AXES order.
            if ~exist('D', 'var'); D = Helmholtz.D_DEFAULT; end
            g = Helmholtz.GAMMA;
            Bq = Bq(:);
            L = zeros(4, 2);
            for k = 1:4
                bz = Helmholtz.NV_AXES(k, :) * Bq;
                bp = sqrt(max(Bq' * Bq - bz^2, 0)) / sqrt(2);
                H = [D + g*bz, g*bp, 0; g*bp, 0, g*bp; 0, g*bp, D - g*bz];
                ev = sort(eig(H));
                L(k, :) = [ev(2) - ev(1), ev(3) - ev(1)];
            end
        end

        function d = dipsAlong(u, rho, D)
            % Where the dips land for |B| = rho along unit direction u (exact Hamiltonian):
            %   d.outer = [lowest highest]                  -- the widest pair
            %   d.inner = [highest lower, lowest upper]     -- the pair nearest D
            % For an NV axis the outer pair is the aligned orientation (D -+ gamma*rho
            % exactly) and the inner pair the three others; along a cube axis all four are
            % equivalent and inner = outer.
            L = Helmholtz.nvLinesByAxis(rho * reshape(u, 1, 3) / norm(u), D);
            d.outer = [min(L(:, 1)), max(L(:, 2))];
            d.inner = [max(L(:, 1)), min(L(:, 2))];
        end

        function rho = rhoForDip(u, D, which, f)
            % The |B| along u that puts dip `which` ('outerLo' 'outerHi' 'innerLo' 'innerHi')
            % at f MHz. Each moves monotonically with |B| below ~120 G, so bisection. NaN if
            % no |B| in 0..120 G puts it there.
            g = @(r) Helmholtz.pickDip(Helmholtz.dipsAlong(u, r, D), which) - f;
            a = 0; b = 120; fa = g(a); fb = g(b);
            rho = NaN;
            if ~isfinite(fa) || ~isfinite(fb) || fa*fb > 0; return; end
            for it = 1:60
                m = 0.5*(a + b); fm = g(m);
                if fa*fm <= 0; b = m; else; a = m; fa = fm; end
            end
            rho = 0.5*(a + b);
        end

        function g = groupDips(f)
            % Distinct dips and how many transitions share each (within 0.5 MHz): rows
            % [frequency multiplicity], by frequency.
            f = sort(f(:));
            g = zeros(0, 2);
            k = 1;
            while k <= numel(f)
                j = k;
                while j < numel(f) && f(j+1) - f(k) < 0.5; j = j + 1; end
                g(end+1, :) = [mean(f(k:j)), j - k + 1];                          %#ok<AGROW>
                k = j + 1;
            end
        end

        function C = plotColours()
            % fitStyle's palette (Setup 7\functions\fitStyle.m), so the GUI reads like the fits.
            C = struct('surface', [252 252 251]/255, 'ink', [20 20 20]/255, ...
                       'ink2', [82 81 78]/255, 's1', [235 104 52]/255, ...
                       's2', [122 61 148]/255, 's3', [107 120 133]/255, ...
                       'ok', [0 115 89]/255, 'band', [219 224 232]/255);
        end

        function styleAxes(ax, C)
            % fitStyle's axes: solid grid at 12% behind the data, ticks out, no box.
            grid(ax, 'on'); box(ax, 'off');
            set(ax, 'FontSize', 9, 'XColor', C.ink2, 'YColor', C.ink2, 'TickDir', 'out', ...
                    'Color', C.surface, 'GridLineStyle', '-', 'GridAlpha', 0.12, ...
                    'Layer', 'bottom');
        end

        function v = pickDip(d, which)
            switch which
                case 'outerLo'; v = d.outer(1);
                case 'outerHi'; v = d.outer(2);
                case 'innerLo'; v = d.inner(1);
                case 'innerHi'; v = d.inner(2);
                otherwise; error('Helmholtz:dip', 'no dip "%s"', which);
            end
        end

        function r = sph2cartPhys(rho, theta, phi)
            % theta POLAR from +z, phi AZIMUTH from +x, degrees. Returns a 1x3 row (or Nx3).
            rho = rho(:); theta = theta(:); phi = phi(:);
            r = [rho .* sind(theta) .* cosd(phi), rho .* sind(theta) .* sind(phi), ...
                 rho .* cosd(theta)];
        end

        function [rho, theta, phi] = cart2sphPhys(r)
            % theta POLAR from +z (0..180), phi AZIMUTH from +x (-180..180), degrees.
            r = reshape(r, [], 3);
            rho = sqrt(sum(r.^2, 2));
            theta = acosd(max(min(r(:, 3) ./ max(rho, realmin), 1), -1));
            theta(rho == 0) = 0;
            phi = atan2d(r(:, 2), r(:, 1));
        end

        function f = nvLines(Bq, D)
            % The 8 transitions (MHz) at field Bq (crystal frame), exact Hamiltonian per
            % orientation: H = D*Sz^2 + gamma*(Bpar*Sz + Bperp*Sx). Sorted.
            if ~exist('D', 'var'); D = Helmholtz.D_DEFAULT; end
            g = Helmholtz.GAMMA;
            Bq = Bq(:);
            f = zeros(1, 8);
            for k = 1:4
                bz = Helmholtz.NV_AXES(k, :) * Bq;
                bp = sqrt(max(Bq' * Bq - bz^2, 0)) / sqrt(2);
                H = [D + g*bz, g*bp, 0; g*bp, 0, g*bp; 0, g*bp, D - g*bz];
                ev = sort(eig(H));
                f(2*k-1) = ev(2) - ev(1);
                f(2*k) = ev(3) - ev(1);
            end
            f = sort(f);
        end

        function [lo, hi] = rayInterval(M, b0, Imin, Imax, u)
            % r such that r*u is reachable, r >= 0. [NaN NaN] if none.
            Minv = M \ eye(3);
            d = (Minv * u(:))';
            e = (-Minv * b0(:))';
            [lo, hi] = Helmholtz.linearInterval(d, e, Imin, Imax);
            lo = max(lo, 0);
            if ~(lo <= hi); lo = NaN; hi = NaN; end
        end

        function [lo, hi] = linearInterval(d, e, Imin, Imax)
            % t such that Imin <= e + t*d <= Imax, every component. The bounds are the exact
            % boundaries; I_TOL only decides whether a coil with no t-dependence is in range.
            tol = Helmholtz.I_TOL;
            lo = -inf; hi = inf;
            for k = 1:3
                if abs(d(k)) < 1e-14
                    if e(k) < Imin(k) - tol || e(k) > Imax(k) + tol
                        lo = NaN; hi = NaN; return
                    end
                    continue
                end
                a = (Imin(k) - e(k)) / d(k);
                b = (Imax(k) - e(k)) / d(k);
                lo = max(lo, min(a, b));
                hi = min(hi, max(a, b));
            end
            if ~(lo <= hi); lo = NaN; hi = NaN; end
        end

        function [lo, hi] = fieldLimits(M, b0, Imin, Imax, Bnow, what)
            what = char(what);   % "Bx" (string) as well as 'Bx'
            % Exact range of one field coordinate with the other two held. M already
            % includes the flip sign. See the class header for the method.
            Minv = M \ eye(3);
            Bnow = reshape(Bnow, 1, 3);
            b0 = reshape(b0, 1, 3);
            switch what
                case {'Bx', 'By', 'Bz'}
                    c = what(2) - 'w';
                    I0 = (Minv * (Bnow - b0)')';
                    d = Minv(:, c)';
                    [tlo, thi] = Helmholtz.linearInterval(d, I0, Imin, Imax);
                    lo = Bnow(c) + tlo; hi = Bnow(c) + thi;
                case 'Brho'
                    [rho, th, ph] = Helmholtz.cart2sphPhys(Bnow);
                    if rho < 1e-12; th = 0; ph = 0; end
                    u = Helmholtz.sph2cartPhys(1, th, ph);
                    [lo, hi] = Helmholtz.rayInterval(M, b0, Imin, Imax, u);
                case {'Btheta', 'Bphi'}
                    [rho, th, ph] = Helmholtz.cart2sphPhys(Bnow);
                    e = (-Minv * b0')';
                    if strcmp(what, 'Btheta')
                        % I = rho*Minv*[sin t cos ph; sin t sin ph; cos t] + e
                        %   = A cos t + C sin t + e
                        A = rho * (Minv * [0; 0; 1])';
                        C = rho * (Minv * [cosd(ph); sind(ph); 0])';
                        [lo, hi] = Helmholtz.angleInterval(A, C, e, Imin, Imax, th, false);
                    else
                        A = rho * sind(th) * (Minv * [1; 0; 0])';
                        C = rho * sind(th) * (Minv * [0; 1; 0])';
                        e = e + rho * cosd(th) * (Minv * [0; 0; 1])';
                        [lo, hi] = Helmholtz.angleInterval(A, C, e, Imin, Imax, ph, true);
                    end
                otherwise
                    error('Helmholtz:limits', 'no coordinate %s', what);
            end
        end

        function [lo, hi] = angleInterval(A, C, e, Imin, Imax, tNow, circular)
            % The stretch of t (degrees) around tNow where Imin <= A cos t + C sin t + e
            % <= Imax holds for all three coils. circular: t is an azimuth (wraps at 360);
            % otherwise t is a polar angle in [0, 180].
            tol = Helmholtz.I_TOL;
            feas = @(t) all(A .* cosd(t) + C .* sind(t) + e >= Imin - tol & ...
                            A .* cosd(t) + C .* sind(t) + e <= Imax + tol);
            % every angle at which some coil crosses a limit: R cos(t - delta) = v
            bp = [];
            for k = 1:3
                R = hypot(A(k), C(k));
                if R < 1e-14; continue; end
                delta = atan2d(C(k), A(k));
                for v = [Imin(k) - e(k), Imax(k) - e(k)]
                    if abs(v) <= R
                        a = acosd(v / R);
                        bp = [bp, delta + a, delta - a];                         %#ok<AGROW>
                    end
                end
            end
            if circular
                bp = unique(round(mod(bp, 360), 10));
                t0 = mod(tNow, 360);
                if isempty(bp)
                    if feas(t0); lo = -180; hi = 180; else; lo = NaN; hi = NaN; end
                    return
                end
                n = numel(bp);
                % segment i runs from bp(i) to bp(i+1) (the last wraps to bp(1) + 360)
                segLo = bp; segHi = [bp(2:end), bp(1) + 360];
                ok = false(1, n);
                for i = 1:n; ok(i) = feas((segLo(i) + segHi(i)) / 2); end
                if all(ok); lo = -180; hi = 180; return; end
                % the feasible segment holding t0 (or touching it, when the present point
                % sits exactly on a supply limit), with t0 unwrapped to lie inside it
                i0 = []; t0u = t0;
                for i = 1:n
                    if ~ok(i); continue; end
                    for cand = [t0, t0 + 360]
                        if cand >= segLo(i) - 1e-9 && cand <= segHi(i) + 1e-9
                            i0 = i; t0u = cand; break
                        end
                    end
                    if ~isempty(i0); break; end
                end
                if isempty(i0); lo = NaN; hi = NaN; return; end
                % grow forwards then backwards through feasible segments, by LENGTH --
                % no modulo arithmetic, so an end that coincides with t0 cannot wrap by 360
                segLen = segHi - segLo;
                total = segLen(i0);
                endU = segHi(i0); j = i0;
                while total < 360 - 1e-9
                    jn = mod(j, n) + 1;
                    if ~ok(jn) || jn == i0; break; end
                    endU = endU + segLen(jn); total = total + segLen(jn); j = jn;
                end
                startU = segLo(i0); j = i0;
                while total < 360 - 1e-9
                    jp = mod(j - 2, n) + 1;
                    if ~ok(jp) || jp == i0; break; end
                    startU = startU - segLen(jp); total = total + segLen(jp); j = jp;
                end
                if total >= 360 - 1e-9; lo = -180; hi = 180; return; end
                lo = tNow + (startU - t0u); hi = tNow + (endU - t0u);
            else
                bp = mod(bp, 360);
                bp = bp(bp > 0 & bp < 180);
                bp = unique(round([0, bp, 180], 10));
                n = numel(bp) - 1;
                ok = false(1, n);
                for i = 1:n; ok(i) = feas((bp(i) + bp(i+1)) / 2); end
                i0 = find(ok & bp(1:n) <= tNow + 1e-9 & bp(2:end) >= tNow - 1e-9, 1);
                if isempty(i0); lo = NaN; hi = NaN; return; end
                iL = i0; while iL > 1 && ok(iL - 1); iL = iL - 1; end
                iR = i0; while iR < n && ok(iR + 1); iR = iR + 1; end
                lo = bp(iL); hi = bp(iR + 1);
            end
        end

        function c = cleanCalibration(c)
            % Shapes as the code expects them, whether from coilMatrixFit or jsondecode
            % (which returns rows of [x,y,z] as columns).
            if isfield(c, 'matrix'); c.matrix = reshape(c.matrix, 3, 3); end
            for f = {'offset', 'offsetErr', 'currentMin', 'currentMax'}
                if isfield(c, f{1}) && numel(c.(f{1})) == 3
                    c.(f{1}) = reshape(c.(f{1}), 1, 3);
                end
            end
            if isfield(c, 'matrixErr') && numel(c.matrixErr) == 9
                c.matrixErr = reshape(c.matrixErr, 3, 3);
            end
            if isfield(c, 'date'); c.date = char(c.date); end
            if isfield(c, 'frame'); c.frame = char(c.frame); end
        end
    end
end

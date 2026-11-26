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

        ESR_WIDTH = 5;          % MHz, half width of the lines in the expected-ESR plot
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
        function plotExpectedESR(obj, gAxes)
            if ~exist('gAxes', 'var')
                figure;
                gAxes = axes(gcf);
            end
            cla(gAxes);
            lines = obj.expectedLines();
            w = obj.ESR_WIDTH;
            span = max(abs(lines - obj.D)) + 6*w;
            f = linspace(obj.D - span, obj.D + span, 1500);
            y = ones(size(f));
            for k = 1:numel(lines)
                y = y - 0.01 * w^2 ./ ((f - lines(k)).^2 + w^2);
            end
            plot(gAxes, f, y, 'k');
            hold(gAxes, 'on');
            for k = 1:numel(lines)
                plot(gAxes, lines(k)*[1 1], [min(y) 1], ':', 'Color', [0.6 0.6 0.6]);
            end
            hold(gAxes, 'off');
            xlim(gAxes, [f(1) f(end)]);
            gAxes.XLabel.String = 'frequency (MHz)';
            gAxes.YLabel.String = 'expected signal';
            if obj.isCalibrated
                gAxes.Title.String = sprintf('|B| = %.2f G, %.2f deg from the nearest NV', ...
                    obj.Brho, acosd(min(obj.nearestNVcos(), 1)));
            else
                gAxes.Title.String = 'UNCALIBRATED: lines assume the coil frame is the crystal frame';
            end
        end

        function plotFieldMap(obj, gAxes)
            % The fields the coils can reach (the parallelepiped), the NV axes, and B.
            if ~exist('gAxes', 'var')
                figure;
                gAxes = axes(gcf);
            end
            cla(gAxes);
            hold(gAxes, 'on');
            % parallelepiped: the 8 corners of the current box mapped through the model
            [a, b, c] = ndgrid([0 1], [0 1], [0 1]);
            frac = [a(:), b(:), c(:)];
            Ic = obj.minCurrent + frac .* (obj.maxCurrent - obj.minCurrent);
            Bc = obj.fieldFor(Ic);
            faces = [1 2 4 3; 5 6 8 7; 1 2 6 5; 3 4 8 7; 1 3 7 5; 2 4 8 6];
            patch(gAxes, 'Vertices', Bc, 'Faces', faces, 'FaceColor', [1 0.6 0.6], ...
                  'FaceAlpha', 0.12, 'EdgeColor', [0.8 0.3 0.3]);
            L = 1.1 * max(sqrt(sum(Bc.^2, 2)));
            dirs = Helmholtz.nvDirections();
            for k = 1:4
                plot3(gAxes, L*[0 dirs(k,1)], L*[0 dirs(k,2)], L*[0 dirs(k,3)], '-', ...
                      'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
                plot3(gAxes, L*[0 dirs(k+4,1)], L*[0 dirs(k+4,2)], L*[0 dirs(k+4,3)], ':', ...
                      'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
                text(gAxes, L*dirs(k,1), L*dirs(k,2), L*dirs(k,3), ...
                     sprintf('[%s]', num2str(round(dirs(k,:)*sqrt(3)))), 'FontSize', 8);
            end
            Bn = obj.B;
            quiver3(gAxes, 0, 0, 0, Bn(1), Bn(2), Bn(3), 0, 'r', 'LineWidth', 2.5, ...
                    'MaxHeadSize', 0.4);
            plot3(gAxes, 0, 0, 0, 'k.', 'MarkerSize', 12);
            hold(gAxes, 'off');
            axis(gAxes, 'equal');
            grid(gAxes, 'on');
            gAxes.XLabel.String = 'B_x (G)';
            gAxes.YLabel.String = 'B_y (G)';
            gAxes.ZLabel.String = 'B_z (G)';
            gAxes.Title.String = obj.frameName;
            view(gAxes, 135, 25);
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

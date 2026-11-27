classdef viewHelmholtz < ViewVBox
    %VIEWHELMHOLTZ The Helmholtz window.
    %
    %   column 1  field control: general, spherical, cartesian, currents
    %   column 2  the expected ESR at this field (and at the one the Align panel is
    %             previewing), interactive
    %   column 3  align B (direction, |B|, dip positions), the coil model in use, save / load
    %   bottom    error messages
    %
    % The diamond-placement controls and drawing are gone: the coil calibration
    % (calibCoilMatrix) measures the crystal's orientation against the coils, so nothing
    % about the diamond's placement is entered by hand any more.

    methods
        % constructor - parent of a main view is the figure itself!
        function obj = viewHelmholtz(parent, controller)
            obj@ViewVBox(parent, controller);

            mainView = ViewHBox(obj, controller);
            column1 = ViewHelmholtzControl(mainView, controller);
            column2 = ViewVBox(mainView, controller);
            column3 = ViewVBox(mainView, controller);
            ViewVBox(mainView, controller);     % dummy column, so column 3 does not end at the window edge

            % column 2: the expected ESR, given the whole column
            viewExpectedESR = ViewExpectedESR(column2, controller);
            column2.setHeights(-1);
            column2width = viewExpectedESR.width;
            column2height = viewExpectedESR.height + 10;

            % column 3
            vNV = ViewNVAlignment(column3, controller);
            vCal = ViewCoilCalibration(column3, controller);
            vSaveLoad = ViewHelmholtzSaveLoad(column3, controller);
            column3.setHeights([vNV.height, vCal.height, vSaveLoad.height]);
            column3width = max([vNV.width, vCal.width, vSaveLoad.width]);
            column3height = vNV.height + vCal.height + vSaveLoad.height + 15;

            dummyWidth = 3;
            columnsWidths = [column1.width, column2width, column3width, dummyWidth];
            mainView.component.Widths = [column1.width, -1, column3width, dummyWidth];
            errorView = ViewError(obj, controller);
            maximumColumnsHeight = max([column1.height, column2height, column3height]);
            obj.component.Heights = [-1, errorView.height];

            obj.width = max([errorView.width, sum(columnsWidths)]) + 20;
            obj.height = errorView.height + maximumColumnsHeight;
        end
    end
end

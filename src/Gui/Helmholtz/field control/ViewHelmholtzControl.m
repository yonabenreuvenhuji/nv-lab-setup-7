classdef ViewHelmholtzControl < ViewVBox
    %VIEWHELMHOLTZCONTROL The left column of the Helmholtz window, top to bottom:
    %   Supplies            control state, output, flip
    %   Field, spherical    |B|, θ, φ
    %   Field, cartesian    Bx, By, Bz
    %   Coil currents       Ix, Iy, Iz
    % Each is a section in the main window's style (MainGuiStyle): its title is a text row,
    % so no title can cover the first row of controls (the clipped "Step:" rows on R2025a+).

    methods
        function obj = ViewHelmholtzControl(parent, controller)
            obj@ViewVBox(parent, controller, 6, MainGuiStyle.GAP);     % (a margin to the window's edge)
            obj.component.BackgroundColor = MainGuiStyle.col().window;

            vGeneral = ViewGeneralControl(obj, controller);
            vBsphere = ViewBsphereControl(obj, controller);
            vBcart = ViewBcartControl(obj, controller);
            vCurrent = ViewCurrentControl(obj, controller);
            uix.Empty('Parent', obj.component);
            boxes = {vGeneral, vBsphere, vBcart, vCurrent};
            heights = cellfun(@(p) p.height, boxes);
            obj.setHeights([heights, -1]);

            obj.height = MainGuiStyle.rowsHeight(heights) + 12;
            obj.width = max(cellfun(@(p) p.width, boxes)) + 12;
        end
    end
end

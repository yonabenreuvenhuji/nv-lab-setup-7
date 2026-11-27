classdef ViewSpcm < GuiComponent
    %VIEWSPCM The detector section of the main window ("Photodiode", or "SPCM" with a photon
    %counter): ViewSpcmCounter in a titled section.

    methods
        function obj = ViewSpcm(parent, controller)
            obj@GuiComponent(parent, controller);
            spcm = getObjByName(Spcm.NAME);
            titleString = 'Photon counter (SPCM)';
            if spcm.hasPhotodiode()
                titleString = 'Photodiode';
            end
            [panel, inner] = MainGuiStyle.section(parent.component, titleString);
            obj.component = panel;
            counterView = ViewSpcmCounter(struct('component', inner), controller, 'isStandalone', false);
            counterView.popoutFcn = @ViewSpcm.popup;
            inner.Heights = -1;
            obj.height = MainGuiStyle.sectionHeight(counterView.height);
            obj.width = counterView.width;
        end
    end

    methods (Static)
        function popup
            GuiControllerSpcmCounter().start;
        end
    end
end

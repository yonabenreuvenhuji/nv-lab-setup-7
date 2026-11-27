classdef ViewPhotodiodeControl < GuiComponent
    %VIEWPHOTODIODECONTROL "Input [Signal (input 1) ▾]": what the photodiode reads, for the
    %live reading and for scans (PhotoDiodeDigitizerNiDaqControlled.currentChannelGUI).
    %
    % Changing it while the reading runs clears the trace (the old points were a different
    % quantity) and re-arms the digitizer for the new inputs.

    properties
        popInput
    end

    properties (Constant)
        ITEMS = {'Signal (input 1)', 'Reference (input 2)', 'Signal & reference', ...
            ['Signal ' char(8722) ' reference'], ['Signal / ' char(10216) 'reference' char(10217)], 'Laser block (input 4)'}
        TIPS = {'The fluorescence photodiode', 'The laser-sampling photodiode', ...
            'Both at once: two readings and their ratio (cannot be used for a scan)', ...
            'Differential: signal minus reference', ...
            'Referenced: signal divided by the reference normalised to its mean (cancels laser noise)', ...
            'The laser-block photodiode'}
    end

    methods
        function obj = ViewPhotodiodeControl(parent, controller)
            obj@GuiComponent(parent, controller);
            box = MainGuiStyle.hbox(parent.component, 4);
            obj.component = box;
            MainGuiStyle.label(box, 'Input', 'left');
            obj.popInput = MainGuiStyle.popup(box, obj.ITEMS, 'Callback', @obj.callbackInput);
            box.Widths = [36 -1];
            obj.height = MainGuiStyle.ROW_H;
            obj.width = 220;
            obj.update()
        end

        function update(obj)
            spcm = getObjByName(Spcm.NAME);
            k = spcm.currentChannelGUI;
            if k >= 1 && k <= numel(obj.ITEMS)
                obj.popInput.Value = k;
                obj.popInput.TooltipString = sprintf('What the photodiode reads, live and in scans.\nNow: %s', obj.TIPS{k});
            end
        end

        function callbackInput(obj, ~, ~)
            spcm = getObjByName(Spcm.NAME);
            newInput = obj.popInput.Value;
            if newInput == spcm.currentChannelGUI; return; end
            spcm.currentChannelGUI = newInput;
            counter = getObjByName(SpcmCounter.NAME);
            if ~isempty(counter)
                if counter.isRunning
                    counter.resetHistory;   % clears the trace, re-arms for the new inputs
                else
                    counter.clearHistory;
                end
            end
            obj.update();
        end
    end
end

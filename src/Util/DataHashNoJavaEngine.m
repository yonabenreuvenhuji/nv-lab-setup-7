classdef DataHashNoJavaEngine < handle
    %DATAHASHNOJAVAENGINE Stand-in for java.security.MessageDigest, for DataHash without Java.
    %
    % From R2026b MATLAB no longer installs a Java runtime, and DataHash's digest engine
    % is java.security.MessageDigest. DataHash keeps using Java when it is available and
    % creates this engine only when it is not. DataHash feeds it exactly the same byte
    % stream (class, size and contents of every element, recursively), so the hash is
    % still content-based: a Sequence edited in place gets a new hash.
    %
    % The digest is two 64-bit keyHash values over that stream (16 bytes), not MD5/SHA,
    % and keyHash values are only guaranteed stable within one MATLAB session. Every
    % DataHash caller in this code uses the result as an in-memory cache key
    % (PulseGenerator's delay-fixed sequence map, the AWG "already loaded" fingerprints),
    % which is all that needs.

    properties (Access = private)
        chunks = {}     % cell array of uint8 row vectors, in update order
    end

    methods
        function update(obj, bytes)
            % Same calling convention as MessageDigest.update(byte[]).
            if ~isa(bytes, 'uint8')
                bytes = typecast(bytes(:), 'uint8');
            end
            obj.chunks{end+1} = reshape(bytes, 1, []);
        end

        function out = digest(obj)
            % Same calling convention as MessageDigest.digest(): returns the hash bytes
            % and resets the engine.
            stream = [obj.chunks{:}];
            if isempty(stream)
                stream = zeros(1, 0, 'uint8');
            end
            h1 = keyHash(stream);
            h2 = keyHash([uint8(255), stream]);     % a second, independent 64 bits
            out = typecast([uint64(h1), uint64(h2)], 'uint8');
            obj.chunks = {};
        end
    end
end

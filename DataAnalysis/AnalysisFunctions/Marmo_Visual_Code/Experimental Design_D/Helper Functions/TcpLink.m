classdef TcpLink < handle
    properties
        sock = -1
        con  = -1
    end
    methods
        function obj = TcpLink(port)
            obj.sock = pnet('tcpsocket', port);
            if obj.sock < 0, error('Could not open port %d', port); end
            pnet(obj.sock, 'setreadtimeout', 1);
            while obj.con < 0                       % wait for Intan to connect
                obj.con = pnet(obj.sock, 'tcplisten');
            end
            pnet(obj.con, 'setreadtimeout', 10);
            pnet(obj.con, 'setwritetimeout', 10);
            [ip, p] = pnet(obj.con, 'gethost');
            fprintf_local = @(varargin) builtin('fprintf', varargin{:});
            fprintf_local('Connected: %s:%d\n', num2str(ip), p);
        end
        function fopen(~), end                      % already open
        function fprintf(obj, fmt, varargin)
            str = sprintf(fmt, varargin{:});
            pnet(obj.con, 'write', [uint8(str) uint8(10)]);   % text + LF
        end
        function s = fscanf(obj, varargin)
            line = pnet(obj.con, 'readline', 1024);
            s = [line newline];                     % helpers do data(1:end-1)
        end
        function flushinput(obj)
            while ~isempty(pnet(obj.con, 'read', 65536, 'uint8', 'noblock')), end
        end
        function fclose(obj)
            pnet(obj.con, 'close'); pnet(obj.sock, 'close');
        end
        function delete(obj)
            try, pnet(obj.con, 'close'); catch, end
            try, pnet(obj.sock, 'close'); catch, end
        end
    end
end
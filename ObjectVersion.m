function value = ObjectVersion(obj)
%value = ObjectVersion(obj) Get version number from class properties
f = fields(obj);
if sum(strcmp('version',f))
     value = obj.version;
elseif sum(strcmp('Version',f))
     value = obj.Version;
elseif sum(strcmp('CurrentVersion',f))
     value = obj.CurrentVersion;
end
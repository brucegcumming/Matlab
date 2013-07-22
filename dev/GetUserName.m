function username = GetUserName()

os = computer;
if strncmp(os,'PCWIN',5)
    username=getenv('username');
elseif strncmp(os,'GLNXA64',7)
    username=getenv('USER');
else
    username=getenv('USER');
end

z = [];
vals = 0:0.01:1;
for x = 0.0:0.01:1
    for y = 0.0:0.01:1
        z = [z x*y];
    end
end
hist(z);
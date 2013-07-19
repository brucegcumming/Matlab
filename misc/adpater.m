l = 1;

pts = [
    0 l*cos(r2d(20)) l*sin(r2d(20));
    0 0 0;
    l*cos(r2d(20))*tan(r2d(25)) l* cos(r2d(20)) l * cos(r2d(20)) * tan(r2d(25))];



acos(dot(pts(1,:), pts(3,:))./(sqrt(sum(pts(3,:).^2)))* sqrt(sum(pts(1,:).^2))) * 180/pi
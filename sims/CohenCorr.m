
%plot correlation vs PD as used in Cohen and Newsome
x = [1:90];
corrs(x) = -0.0011*x+0.22;
x = [91:135];
corrs(x) = -0.0011*x+0.19;
x = [135:180];
corrs(x) = 0.00156*x -0.17;
plot([1:180],corrs);
hold on;
%superimpose approx values from Fig 4 of Cohen and Newsome 08
plot([22.5 67.5 112.5 157.5],[0.2 0.16 0.12 0.04; 0.16 0.1 0.07 0.12],'o');
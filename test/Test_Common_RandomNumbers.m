% Check the correctness of the random numbers output from Test_Common.ex

% UniformRandom01 (check qqplot)
filename = 'Test_Common_UniformRandom01.txt';
a = load(filename);
b = random('Uniform', 0, 1, size(a));
figure(1);
qqplot(a, b); title('Uniform');
% Note: The following code does not work. There is something wrong
% in matlab. pd has no field 'DistName' needed by qqplot
%pd = makedist('Uniform');
%qqplot(a, pd);

% StandardNormal (check qqplot)
filename = 'Test_Common_StandardNormal.txt';
a = load(filename);
b = random('Normal', 0, 1, size(a));
figure(2);
qqplot(a, b); title('Normal');

% StudentT1 (check qqplot)
filename = 'Test_Common_StudentT1.txt';
a = load(filename);
b = random('T', 1, size(a));
tail = ceil(length(a)/100); % Exclude the tail
a = sort(a); a = a(tail+1:end-tail);
b = sort(b); b = b(tail+1:end-tail);
figure(3);
qqplot(a, b); title('T (df = 1)');

% MultivariateStudentT1 (check marginal)
filename = 'Test_Common_MultivariateStudentT1.txt';
a = load(filename);
b = random('T', 1, [size(a,1),1]);
tail = ceil(size(a,1)/100); % Exclude the tail
a1 = a(:,1);
a2 = a(:,2);
a1 = sort(a1); a1 = a1(tail+1:end-tail);
a2 = sort(a2); a2 = a2(tail+1:end-tail);
b = sort(b); b = b(tail+1:end-tail);
figure(4);
qqplot(a1, b); title('Multivariate T (df = 1), dim 1');
figure(5);
qqplot(a2, b); title('Multivariate T (df = 1), dim 2');

% RandomSech (check pdf)
filename = 'Test_Common_RandomSech.txt';
a = load(filename);
[f, xi] = ksdensity(a);
figure(6);
plot(xi, f); title('KDE, hyperpolic secant');

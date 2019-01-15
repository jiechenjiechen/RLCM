% Plot the relative residual results output from Test_PCG.ex

relres1 = load('Test_PCG.relres1');
relres2 = load('Test_PCG.relres2');
relres3 = load('Test_PCG.relres3');
relres4 = load('Test_PCG.relres4');
semilogy(relres1, 'r');
hold on;
semilogy(relres2, 'b');
semilogy(relres3, 'g');
semilogy(relres4, 'm');
hold off;
legend('No preconditioner', ...
       'No precond, actual res', ...
       'With preconditioner', ...
       'With precond, actual res', ...
       'location', 'southwest');


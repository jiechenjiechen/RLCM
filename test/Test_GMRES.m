% Plot the relative residual results output from Test_GMRES.ex

relres1 = load('Test_GMRES.relres1');
relres2 = load('Test_GMRES.relres2');
relres3 = load('Test_GMRES.relres3');
relres4 = load('Test_GMRES.relres4');
semilogy(relres1, 'r');
hold on;
semilogy(relres2, 'b');
semilogy(relres3, 'g');
semilogy(relres4, 'm');
hold off;
legend('No preconditioner, restart = 100', ...
       'With preconditioner, restart = 100', ...
       'No preconditioner, restart = 10', ...
       'With preconditioner, restart = 10', ...
       'location', 'southwest');


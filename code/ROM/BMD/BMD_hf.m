% BMD model of size r generated here

function[F_mat, G_mat, H_mat, D_mat, L_mat, E_mat, W_mat, V, rel_fit_error_vec] = BMD_hf(r, X_0, X_1, U_0, Y_0, U_1, L_c_mat, L_o_mat)
[n_x, n_s, n_g] = size(X_0);
[n_u,~,~] = size(U_0);
[n_y,~,~] = size(Y_0);

% Building projection-matrices W_mat
Q_bar = zeros(n_x, n_g*r);
for k = 1:n_g
    L_c = L_c_mat{k};
    L_o = L_o_mat{k};
    [U, Sigma, ~] = svd(L_c.'*L_o);
    [U_bar,~,~] = svd(L_c*U(:,1:r));
    Q_bar(:,1 + r*(k-1):k*r) = U_bar(:,1:r);
end

[Q,~,~] = svd(Q_bar);
V = Q(:,1:r); % V (basis space, or primal modes in balanced truncation literature) is constant over all velocities
W_mat = zeros(n_x, r, n_g); % W (test space, or adjoint modes in balanced truncation literature) is velocity dependent (parameter varying)
for k = 1:n_g
    [Q, R] = qr(L_o_mat{k}.'*V);
    Q1 = Q(:,1:r); % Thin QR
    R1 = R(1:r,:);
    W_mat(:,:,k) = L_o_mat{k}*Q1*(R1.')^-1;
end

% Calculating the local models
F_mat = zeros(r, r, n_g);
G_mat = zeros(r, n_u, n_g);
H_mat = zeros(n_y, r, n_g);
D_mat = zeros(n_y, n_u, n_g);
L_mat = zeros(r, n_u, n_g);
E_mat = zeros(n_y, n_u, n_g);
rel_fit_error_vec = zeros(1, n_g);
for k = 1:n_g    
    ss_matrices = [W_mat(:,:,k).'*X_1(:,:,k); Y_0(:,:,k)]*pinv([W_mat(:,:,k).'*X_0(:,:,k); U_0(:,:,k); U_1(:,:,k)]);
    F_mat(:,:,k) = ss_matrices(1:r, 1:r);
    G_mat(:,:,k) = ss_matrices(1:r, (r+1):(r+1+n_u-1));
    L_mat(:,:,k) = ss_matrices(1:r, (r+1+n_u):end);
    H_mat(:,:,k) = ss_matrices((r+1):end, 1:r);
    D_mat(:,:,k) = ss_matrices((r+1):end, (r+1):(r+1+n_u-1));
    E_mat(:,:,k) = ss_matrices((r+1):end, (r+1+n_u):end);
    
    % Check the shapes of the matrices (unnecessary)
    assert(size(F_mat,1) == r && size(F_mat,2) == r);
    assert(size(G_mat,1) == r && size(G_mat,2) == n_u);
    assert(size(H_mat,1) == n_y && size(H_mat,2) == r);
    assert(size(D_mat,1) == n_y && size(D_mat,2) == n_u);
    assert(size(L_mat,1) == r && size(L_mat,2) == n_u);
    assert(size(E_mat,1) == n_y && size(E_mat,2) == n_u);
    
    e_term_1 = [W_mat(:,:,k).'*X_1(:,:,k); Y_0(:,:,k)] - ss_matrices*[W_mat(:,:,k).'*X_0(:,:,k); U_0(:,:,k); U_1(:,:,k)];
    e_term_2 = X_1(:,:,k);
    e_term_3 = W_mat(:,:,k).'*X_1(:,:,k);
    rel_fit_error_vec(k) = norm(e_term_1, 'fro').^2 + norm(e_term_2, 'fro').^2 ...
        - norm(e_term_3, 'fro').^2;
    rel_fit_error_vec(k) = rel_fit_error_vec(k)./norm([X_1(:,:,k); Y_0(:,:,k)], 'fro');
end

end
%% my parameter-varying aIODMD function giving back the matrices
% r: rank for truncation (if r = 0, we perform optimal way)
% X_0: state matrix: n_x x Nt x n_g
% X_1: shifted state matrix: n_x x Nt x n_g
% U_0: input matrix: n_u x n_g
% Y_0: output matrix: n_y x n_g
function[F_mat, G_mat, H_mat, D_mat, L_mat, E_mat, Q, rel_fit_error_vec] = aIODMD_hf(r, X_0, X_1, U_0, Y_0, U_1)
[n_x, n_s, n_g] = size(X_0);
[n_u,~,~] = size(U_0);
[n_y,~,~] = size(Y_0);

% reshape to column
X_0_re = reshape(X_0, [n_x, n_s*n_g]);
X_1_re = reshape(X_1, [n_x, n_s*n_g]);
U_0_re = reshape(U_0, [n_u, n_s*n_g]);
Y_0_re = reshape(Y_0, [n_y, n_s*n_g]);

% Choose the global Q-matrix
[U, Sigma_glob, V_glob] = svd(X_0_re, 'econ'); % reshaped data-matrix

U_r = U(:,1:r);
Q = U_r;

% Calculating the local models
F_mat = zeros(r, r, n_g);
G_mat = zeros(r, n_u, n_g);
H_mat = zeros(n_y, r, n_g);
D_mat = zeros(n_y, n_u, n_g);
L_mat = zeros(r, n_u, n_g);
E_mat = zeros(n_y, n_u, n_g);
rel_fit_error_vec = zeros(1, n_g);
for k = 1:n_g
    ss_matrices = [Q.'*X_1(:,:,k); Y_0(:,:,k)]*pinv([Q.'*X_0(:,:,k); U_0(:,:,k); U_1(:,:,k)]);
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
    
    e_term_1 = [Q.'*X_1(:,:,k); Y_0(:,:,k)] - ss_matrices*[Q.'*X_0(:,:,k); U_0(:,:,k); U_1(:,:,k)];
    e_term_2 = X_1(:,:,k);
    e_term_3 = Q.'*X_1(:,:,k);
    rel_fit_error_vec(k) = norm(e_term_1, 'fro').^2 + norm(e_term_2, 'fro').^2 ...
        - norm(e_term_3, 'fro').^2;
    rel_fit_error_vec(k) = rel_fit_error_vec(k)./norm([X_1(:,:,k); Y_0(:,:,k)], 'fro');
end

end
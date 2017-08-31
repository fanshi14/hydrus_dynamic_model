syms link_weight_vec sum_weight n_links link_length D R_local T_local
syms link_weight_1 link_weight_2 link_weight_3 link_weight_4
syms px py pz er ep ey q1 q2 q3
syms D11 D12 D13 D22 D23 D33
D = sym(zeros(9, 9));
D11 = sym(zeros(3, 3));
D12 = sym(zeros(3, 3));
D13 = sym(zeros(3, 3));
D22 = sym(zeros(3, 3));
D23 = sym(zeros(3, 3));
D33 = sym(zeros(3, 3));
link_weight_vec = [link_weight_1; link_weight_2; link_weight_3; link_weight_4];
sum_weight = link_weight_1 + link_weight_2 + link_weight_3 + link_weight_4;
syms D_origin
D_origin = sym(zeros(9, 9));
%% debug
load_mid_result_flag = true;
load_simplify_C_mid_result_flag = true;
if load_mid_result_flag
    disp('Load mid reuslt.')
else
    disp('Do NOT load mid reuslt.')
end
disp('start time:')
disp(datestr(now))

syms H_inertial
H_inertial = sym(zeros(3, 3, 4));
for i = 1:4
    H_inertial(1, 1, i) = 0.0001;
    H_inertial(2, 2, i) = 0.0001;
    H_inertial(3, 3, i) = 0.0002;
end

R_local = sym(zeros(3, 3));
T_local = sym(zeros(3, 3));
R_local = [cos(ey) -sin(ey) 0; sin(ey) cos(ey) 0; 0 0 1] * ...
          [cos(ep) 0 sin(ep); 0 1 0; -sin(ep) 0 cos(ey)] * ...
          [1 0 0; 0 cos(er) -sin(er); 0 sin(er) cos(er)];
%% convert euler rate to angular velocity (base frame)
T_local = [1 0 -sin(ep); 0 cos(er) cos(ep)*sin(er); 0 -sin(er) ...
           cos(ep)*cos(er)];

syms S_operation(x, y, z) %% skew_symmetric matrix
S_operation(x, y, z) = [0 -z y; z 0 -x; -y x 0 ];

%% manipulator
syms link_center_pos_local_vec link_end_pos_local_vec
link_center_pos_local_vec = sym(zeros(3, 4));
link_end_pos_local_vec = sym(zeros(3, 4));

link_center_pos_local_vec(1, 1) = link_length / 2.0;
link_end_pos_local_vec(1, 1) = link_length;
%% todo: currently only consider one dof in every joint
syms R_li_b
R_li_b = sym(zeros(3, 3, 4));
R_li_b(:, :, 1) = sym(eye(3));
R_li_b(:, :, 2) = [cos(q1) -sin(q1) 0; sin(q1) cos(q1) 0; 0 0 1];
R_li_b(:, :, 3) = [cos(q1+q2) -sin(q1+q2) 0; sin(q1+q2) cos(q1+q2) 0; 0 0 1];
R_li_b(:, :, 4) = [cos(q1+q2+q3) -sin(q1+q2+q3) 0; sin(q1+q2+q3) cos(q1+q2+q3) ...
        0; 0 0 1];
%% FK
for i = 2:4
    link_center_pos_local_vec(:, i) = link_end_pos_local_vec(:, i-1) + ...
        R_li_b(:, :, i) * [link_length/2.0; 0; 0];
    link_end_pos_local_vec(:, i) = link_end_pos_local_vec(:, i-1) + ...
        R_li_b(:, :, i) * [link_length; 0; 0];
end
%% Jacobian
syms Jacobian_P Jacobian_W
Jacobian_P = sym(zeros(3, 4, 4));
Jacobian_W = sym(zeros(3, 4, 4));
for i = 1:4
    %% todo: currently only consider one dof in every joint
    z_axis = sym([0; 0; 1]);
    dist_to_end = sym(zeros(3, 1));
    for j = 1:i
        if j == 1
            dist_to_end = link_center_pos_local_vec(:, i);
        else
            dist_to_end = link_center_pos_local_vec(:, i) - ...
                link_end_pos_local_vec(:, j-1);
        end
        Jacobian_P(:, j, i) = cross(z_axis, dist_to_end);
        Jacobian_W(:, j, i) = z_axis;
    end
end

Jacobian_P = simplify(Jacobian_P); %% simplify matrix
%% There is no joint before root link
Jacobian_P = Jacobian_P(:, 2:4, :);
Jacobian_W = Jacobian_W(:, 2:4, :);

%% calculate D
if not(load_mid_result_flag)
    for i = 1:3
        D11(i, i) = sum_weight;
    end
    D_origin(1:3, 1:3) = D11;
    D(1:3, 1:3) = D11;

    for i = 1:4
        p_bli = R_local * link_center_pos_local_vec(:, i);
        S_p_bli = S_operation(p_bli(1), p_bli(2), p_bli(3));
        D12 = D12 - link_weight_vec(i) * S_p_bli * T_local;
    end
    D_origin(1:3, 4:6) = D12;
    D_origin(4:6, 1:3) = (D12.');
    D12 = simplify(D12); %% simplify matrix
    D(1:3, 4:6) = D12;
    D(4:6, 1:3) = (D12.');

    for i = 1:4
        D13 = D13 + link_weight_vec(i) * R_local * Jacobian_P(:, :, i);
    end
    D_origin(1:3, 7:9) = D13;
    D_origin(7:9, 1:3) = (D13.');
    D13 = simplify(D13); %% simplify matrix
    D(1:3, 7:9) = D13;
    D(7:9, 1:3) = (D13.');

    for i = 1:4
        p_bli = R_local * link_center_pos_local_vec(:, i);
        S_p_bli = S_operation(p_bli(1), p_bli(2), p_bli(3));
        D22 = D22 + link_weight_vec(i) * (T_local.') * (S_p_bli.') * ...
              S_p_bli * T_local + (T_local.') * R_local * R_li_b(:, :, i) * ...
              H_inertial(:, :, i) * (R_li_b(:, :, i).') * (R_local.') * T_local;
    end
    D_origin(4:6, 4:6) = D22;
    D22 = simplify(D22); %% simplify matrix
    D(4:6, 4:6) = D22;

    for i = 1:4
        p_bli = R_local * link_center_pos_local_vec(:, i);
        S_p_bli = S_operation(p_bli(1), p_bli(2), p_bli(3));
        D23 = D23 + (T_local.') * R_local * R_li_b(:, :, i) * H_inertial(:, :, i) * ...
              (R_li_b(:, :, i).') * Jacobian_W(:, :, i)...
              - link_weight_vec(i) * (T_local.') * (S_p_bli.') * R_local * ...
              Jacobian_P(:, :, i);
    end
    D_origin(4:6, 7:9) = D23;
    D_origin(7:9, 4:6) = (D23.');
    D23 = simplify(D23); %% simplify matrix
    D(4:6, 7:9) = D23;
    D(7:9, 4:6) = (D23.');

    for i = 1:4
        D33 = D33 + link_weight_vec(i) * (Jacobian_P(:, :, i).') * Jacobian_P(:, :, i) ...
              + (Jacobian_W(:, :, i).') * R_li_b(:, :, i) * H_inertial(:, :, i) ...
              * (R_li_b(:, :, i).') * Jacobian_W(:, :, i);
    end
    D_origin(7:9, 7:9) = D33;
    D33 = simplify(D33); %% simplify matrix
    D(7:9, 7:9) = D33;

    disp('D is generated.');
end

syms q_vec d_q_vec
q_vec = [px; py; pz; er; ep; ey; q1; q2; q3];
syms d_px d_py d_pz d_er d_ep d_ey d_q1 d_q2 d_q3
d_q_vec = [d_px; d_py; d_pz; d_er; d_ep; d_ey; d_q1; d_q2; d_q3];
syms C
C = sym(zeros(9, 9));
%% calculate C
if not(load_mid_result_flag)
    for k = 1:9
        for j = 1:9
            for i = 1:9
                C(k, j) = C(k, j) + (diff(D(k, j), q_vec(i)) + ...
                                     diff(D(k, i), q_vec(j)) - ...
                                     diff(D(i, j), q_vec(k))) ...
                          * 0.5 * d_q_vec(i);
            end
        end
    end
    %C = simplify(C); %% simplify matrix
    disp('C is generated.');
end

syms g
%% calculate g
if not(load_mid_result_flag)
    syms U
    U = 0;
    for i = 1:4
        U = U + link_weight_vec(i) * 9.78 *  ([0;0;1].') ...
            * ([px; py; pz] + R_local * link_center_pos_local_vec(:, i));
    end
    U = simplify(U); %% simplify matrix

    g = sym(zeros(9, 1));
    for i = 1:9
        g(i) = diff(U, q_vec(i));
    end
    g = simplify(g); %% simplify matrix
    disp('g is generated.');
end

syms B B1 B2
syms f1 f2 f3 f4 tau1 tau2 tau3
syms u
u = [f1; f2; f3; f4; tau1; tau2; tau3];
syms P_cog_local
P_cog_local = sym(zeros(3, 1));
for i = 1:4
    P_cog_local = P_cog_local + link_weight_vec(i) * ...
        link_center_pos_local_vec(:, i);
end
P_cog_local = P_cog_local / sum_weight;

%% calculate B (B is the funciton of B(u))
if not(load_mid_result_flag)
    B = sym(zeros(9, 1));
    B1 = sym(zeros(9, 9));
    B2 = sym(zeros(9, 1));
    B1(1:3, 1:3) = R_local;
    B1(4:6, 4:6) = (T_local.') * R_local;
    B1(7:9, 7:9) = sym(eye(3));
    B2(7:9, 1) = u(5:7, 1);

    syms momentum_local
    momentum_local = sym(zeros(3, 1));
    for i = 1:4
        momentum_local = momentum_local + cross(link_center_pos_local_vec(:, i), [0; 0; u(i)]);
    end

    B2(4:6, 1) = momentum_local;
    B2(1:3, 1) = R_local * [0; 0; f1 + f2 + f3 + f4];
    B = B1 * B2;
    B = simplify(B); %% simplify matrix
    disp('B is generated.');
end

%% load mid result
if load_mid_result_flag
    load('hydrus_mid_result.mat');
    disp('D, C, g, B is loaded.')
%% save mid result
else
    save('hydrus_mid_result.mat', 'D', 'D_origin', 'C', 'g', 'B');
    disp('D, C, g, B data is saved.');
end

%% load simplify C
syms C1
if load_simplify_C_mid_result_flag
    load('hydrus_mid_result_simple_C.mat');
    C = C1;
    disp('load simplify C mid result.');
end

%% D * s'' + C * s' + g = B

%% simplified state and control input
syms qs_vec us_vec
qs_vec = [px; py; pz; er; ep; ey; ...
          d_px; d_py; d_pz; d_er; d_ep; d_ey];
us_vec = [f1; f2; f3; f4];
syms Ds Cs gs Bs
Ds = D(1:6, 1:6);
Cs = C(1:6, 1:6);
gs = g(1:6, 1);
Bs = B(1:6, 1);

syms Csds Ds3 Cs3
Csds = Cs * qs_vec(7:12, 1);
Ds3 = D(1:6, 7:9);
Cs3 = C(1:6, 7:9);
syms Ds_x Bs_x Bs_u Csds_x Csds_dx gs_x Ds3_x Cs3_x Cs3_dx
Ds_x = sym(zeros(6, 6, 6));
Bs_x = sym(zeros(6, 1, 6));
Bs_u = sym(zeros(6, 1, 4));
Csds_x = sym(zeros(6, 1, 6));
Csds_dx = sym(zeros(6, 1, 6));
gs_x = sym(zeros(6, 1, 6));
Ds3_x = sym(zeros(6, 3, 6));
Cs3_x = sym(zeros(6, 3, 6));
Cs3_dx = sym(zeros(6, 3, 6));

for i = 1:6
    Ds_x(:, :, i) = diff(Ds, qs_vec(i));
    Bs_x(:, :, i) = diff(Bs, qs_vec(i));
    Csds_x(:, :, i) = diff(Csds, qs_vec(i));
    Csds_dx(:, :, i) = diff(Csds, qs_vec(i+6));
    gs_x(:, :, i) = diff(gs, qs_vec(i));
    Ds3_x(:, :, i) = diff(Ds3, qs_vec(i));
    Cs3_x(:, :, i) = diff(Cs3, qs_vec(i));
    Cs3_dx(:, :, i) = diff(Cs3, qs_vec(i+6));
end

for i = 1:4
    Bs_u(:, :, i) = diff(Bs, us_vec(i));
end



%% save Ds Ds3 C Cs3 Csds gs Bs ...
%%      Ds_x Bs_x Bs_u Csds_x Csds_dx gs_x Ds3_x Cs3_x Cs3_dx
fid = fopen('hrdrus_matrix_Ds.txt', 'wt');
var = Ds; %Ds
[row, col] = size(var);
fprintf(fid, '\n\n\n Ds, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('Ds output is finished.');

fid = fopen('hrdrus_matrix_Ds3.txt', 'wt');
var = Ds3; %Ds3
[row, col] = size(var);
fprintf(fid, '\n\n\n Ds3, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('Ds3 output is finished.');

fid = fopen('hrdrus_matrix_C.txt', 'wt');
var = C; %C
[row, col] = size(var);
fprintf(fid, '\n\n\n C, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('C output is finished.');

fid = fopen('hrdrus_matrix_Cs3.txt', 'wt');
var = Cs3; %Cs3
[row, col] = size(var);
fprintf(fid, '\n\n\n Cs3, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('Cs3 output is finished.');

fid = fopen('hrdrus_matrix_Csds.txt', 'wt');
var = Csds; %Csds
[row, col] = size(var);
fprintf(fid, '\n\n\n Csds, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('Csds output is finished.');

fid = fopen('hrdrus_matrix_gs.txt', 'wt');
var = gs; %gs
[row, col] = size(var);
fprintf(fid, '\n\n\n gs, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('gs output is finished.');

fid = fopen('hrdrus_matrix_Bs.txt', 'wt');
var = Bs; %Bs
[row, col] = size(var);
fprintf(fid, '\n\n\n Bs, size: %d * %d\n', row, col);
for i = 1:row
    for j = 1:col
        fprintf(fid, '%s\n\n', ccode(var(i,j)));
    end
    fprintf(fid, '\n\n');
end
fclose(fid);
disp('Bs output is finished.');


fid = fopen('hrdrus_matrix_Ds_x.txt', 'wt');
var = Ds_x; %Ds_x
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Ds_x:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Ds_x output is finished.');

fid = fopen('hrdrus_matrix_Bs_x.txt', 'wt');
var = Bs_x; %Bs_x
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Bs_x:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Bs_x output is finished.');

fid = fopen('hrdrus_matrix_Bs_u.txt', 'wt');
var = Bs_u; %Bs_u
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Bs_u:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Bs_u output is finished.');

fid = fopen('hrdrus_matrix_Csds_x.txt', 'wt');
var = Csds_x; %Csds_x
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Csds_x:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Csds_x output is finished.');

fid = fopen('hrdrus_matrix_Csds_dx.txt', 'wt');
var = Csds_dx; %Csds_dx
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Csds_dx:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Csds_dx output is finished.');

fid = fopen('hrdrus_matrix_gs_x.txt', 'wt');
var = gs_x; %gs_x
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n gs_x:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('gs_x output is finished.');

fid = fopen('hrdrus_matrix_Ds3_x.txt', 'wt');
var = Ds3_x; %Ds3_x
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Ds3_x:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Ds3_x output is finished.');

fid = fopen('hrdrus_matrix_Cs3_x.txt', 'wt');
var = Cs3_x; %Cs3_x
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Cs3_x:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Cs3_x output is finished.');

fid = fopen('hrdrus_matrix_Cs3_dx.txt', 'wt');
var = Cs3_dx; %Cs3_dx
[row, col, num] = size(var);
for k = 1:num
    fprintf(fid, '\n\n\n Cs3_dx:[%d], size: %d * %d\n', num, row, col);
    for i = 1:row
        for j = 1:col
            fprintf(fid, '%s\n\n', ccode(var(i,j, num)));
        end
        fprintf(fid, '\n\n');
    end
end
fclose(fid);
disp('Cs3_dx output is finished.');


disp('finish time:')
disp(datestr(now))
disp('Finished.');

syms link_weight_vec sum_weight n_links link_length D R_local T_local
syms link_weight_1 link_weight_2 link_weight_3 link_weight_4
syms px py pz er ep ey q0 q1 q2
syms D11 D12 D13 D22 D23 D33
D = sym(zeros(10, 10));
D11 = sym(zeros(3, 3));
D12 = sym(zeros(3, 3));
D13 = sym(zeros(3, 4));
D22 = sym(zeros(3, 3));
D23 = sym(zeros(3, 4));
D33 = sym(zeros(4, 4));
link_weight_vec = [link_weight_1; link_weight_2; link_weight_3; link_weight_4];

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
R_li_b(:, :, 2) = [cos(q0) -sin(q0) 0; sin(q0) cos(q0) 0; 0 0 1];
R_li_b(:, :, 3) = [cos(q0+q1) -sin(q0+q1) 0; sin(q0+q1) cos(q0+q1) 0; 0 0 1];
R_li_b(:, :, 4) = [cos(q0+q1+q2) -sin(q0+q1+q2) 0; sin(q0+q1+q2) cos(q0+q1+q2) ...
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
            dist_to_end = link_center_pos_local_vec(:, 1);
        else
            dist_to_end = link_center_pos_local_vec(:, i) - ...
                link_end_pos_local_vec(:, j-1);
        end
        Jacobian_P(:, j, i) = cross(z_axis, dist_to_end);
        Jacobian_W(:, j, i) = z_axis;
    end
end

for i = 1:3
    D11(i, i) = sum_weight;
end

for i = 1:4
    p_bli = R_local * link_center_pos_local_vec(:, i);
    S_p_bli = S_operation(p_bli(1), p_bli(2), p_bli(3));
    D12 = D12 - link_weight_vec(i) * S_p_bli * T_local;
end

for i = 1:4
    D13 = D13 + link_weight_vec(i) * R_local * Jacobian_P(:, :, i);
end

for i = 1:4
    p_bli = R_local * link_center_pos_local_vec(:, i);
    S_p_bli = S_operation(p_bli(1), p_bli(2), p_bli(3));
    D22 = D22 + link_weight_vec(i) * (T_local.') * (S_p_bli.') * ...
          S_p_bli * T_local + (T_local.') * R_local * R_li_b(:, :, i) * ...
          H_inertial(:, :, i) * (R_li_b(:, :, i).') * (R_local.') * T_local;
end

for i = 1:4
    p_bli = R_local * link_center_pos_local_vec(:, i);
    S_p_bli = S_operation(p_bli(1), p_bli(2), p_bli(3));
    D23 = D23 + (T_local.') * R_local * R_li_b(:, :, i) * H_inertial(:, :, i) * ...
          (R_li_b(:, :, i).') * Jacobian_W(:, :, i)...
    - link_weight_vec(i) * (T_local.') * (S_p_bli.') * R_local * ...
        Jacobian_W(:, :, i);
end

for i = 1:4
    D33 = D33 + link_weight_vec(i) * (Jacobian_P(:, :, i).') * Jacobian_P(:, :, i) ...
          + (Jacobian_W(:, :, i).') * R_li_b(:, :, i) * H_inertial(:, :, i) ...
          * (R_li_b(:, :, i).') * Jacobian_W(:, :, i);
end

D(1:3, 1:3) = D11;
D(1:3, 4:6) = D12;
D(4:6, 1:3) = (D12.');
D(1:3, 7:10) = D13;
D(7:10, 1:3) = (D13.');
D(4:6, 4:6) = D22;
D(4:6, 7:10) = D23;
D(7:10, 4:6) = (D23.');
D(7:10, 7:10) = D33;
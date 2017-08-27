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
rot1 = [cos(q0) -sin(q0) 0; sin(q0) cos(q0) 0; 0 0 1];
rot2 = [cos(q0+q1) -sin(q0+q1) 0; sin(q0+q1) cos(q0+q1) 0; 0 0 1];
rot3 = [cos(q0+q1+q2) -sin(q0+q1+q2) 0; sin(q0+q1+q2) cos(q0+q1+q2) ...
        0; 0 0 1];
%% FK
link_center_pos_local_vec(:, 2) = link_end_pos_local_vec(:, 1) + ...
    rot1 * [link_length/2.0; 0; 0];
link_end_pos_local_vec(:, 2) = link_end_pos_local_vec(:, 1) + ...
    rot1 * [link_length; 0; 0];
link_center_pos_local_vec(:, 3) = link_end_pos_local_vec(:, 2) + ...
    rot2 * [link_length/2.0; 0; 0];
link_end_pos_local_vec(:, 3) = link_end_pos_local_vec(:, 2) + ...
    rot2 *[link_length; 0; 0];
link_center_pos_local_vec(:, 4) = link_end_pos_local_vec(:, 3) + ...
    rot3 * [link_length/2.0; 0; 0];
link_end_pos_local_vec(:, 4) = link_end_pos_local_vec(:, 3) + ...
    rot3 * [link_length; 0; 0];
%% Jacobian
syms Jacobian_P
Jacobian_P = sym(zeros(3, 4, 4));
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
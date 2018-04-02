function q = quaternions_correctsign(q, q_predicted)

% now we check if we want q or -q
plus  = q_predicted - q;
minus = q_predicted + q;
if (norm(plus) > norm(minus))
    q = - q;
end

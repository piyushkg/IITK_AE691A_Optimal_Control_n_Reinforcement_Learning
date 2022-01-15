function av_act =  possible_actions(state, R)
    curr_state_row = squeeze(R(state(1),state(2),:));
    av_act = find(curr_state_row > -1);
end
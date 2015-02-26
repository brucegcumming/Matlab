function [E,V, pc] = SpikePCA(spk, probes, ids)

if length(probes) == 2
spks = [spk{probes(1)}.values(ids{1},:) spk{probes(2)}.values(ids{2},:)];
[E,V] = eig(cov(spks));
pc(1,:) = (E(:,end)'-mean(E(:,end))) * spks';
pc(2,:) = (E(:,end-1)'-mean(E(:,end-1))) * spks';
pc(3,:) = (E(1:32,end)'-mean(E(1:32,end-1))) * spks(:,1:32)';
pc(4,:) = (E(33:end,end)'-mean(E(33:end,end-1))) * spks(:,33:end)';
pc(5,:) = (E(1:32,end-1)'-mean(E(1:32,end-1))) * spks(:,1:32)';
pc(6,:) = (E(33:end,end-1)'-mean(E(33:end,end-1))) * spks(:,33:end)';
end




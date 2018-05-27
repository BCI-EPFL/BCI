
function xAfter = change_structure_period(xBefore)

xAfter = [];
for iPeriod = 1:size(xBefore,1)
    t = squeeze(xBefore(iPeriod,:,:));
    xAfter= [xAfter ; t];
end

end

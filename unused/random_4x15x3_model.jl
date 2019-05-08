# 4x5x3 model

wet3d = falses(4,15,3)
wet3d[2, :, 1] .= true 
wet3d[2, :, 2] .= true 
wet3d[3:4, 7:9, 1] = true 
wet3d[3:4, 12:14, 1] = true 
wet3d[3, 7:9, 2] = true 
wet3d[3, 12:14, 2] = true 
wet3d[2, 2:4, 3] = true 
wet3d[2, 7:9, 3] = true 
wet3d[3, 12:14, 3] = true 
wet3d[3, 7, 3] = true 
wet3d[3, 9, 3] = true 
wet3d[2, 12, 3] = true 

const iwet = (LinearIndices(wet3d))[findall(!iszero, wet3d)] # replaces find(wet3d) :(nwet = length(iwet)




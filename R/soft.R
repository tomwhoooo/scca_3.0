soft <-
function(y,thr){
    sign(y)*(abs(y)-thr)*(abs(y)>thr)
}


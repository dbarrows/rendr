library(db)
load_all()

sys <- rsys(
        network = network('A + B -> 0, 2e-1'),
        state = c(100, 50),
        T = 0.5
    )
L <- 5
dsys <- rdsys(
        network = sys$network,
        volume = {
            v <- volume(dims = c(L, L, 1), h = 1, c(0, 0))
            #volume_set(v, c(ceiling(L/2), ceiling(L/2), 1), sys$state)
            volume_set(v, c(1, 1, 1), sys$state)
            v
        },
        D = rep(L/sys$T, 2),
        T = 0.5
    )
sol <- dsys |> nsm()
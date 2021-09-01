library(db)
load_all()

dsys <- rdsys(
        network = network('
            A -> 2A, 1
            A + B -> 2B, 5e-3
            B -> 0, 6e-1
        '),
        volume = {
            v <- volume(dims = c(2, 1, 1), h = 1, c(0, 0))
            volume_set(v, c(1, 1, 1), c(50, 100))
            v
        },
        D = c(1, 1),
        T = 10
    )

trajectories <- 1e4

## standard averaged trajectories
# sol1 <- dsys |>
#     nsm(trajectories = trajectories, average = TRUE, parallel = TRUE)
# sol1$u[[100]]

sol2 <- nsm_cpp_pest(compile(dsys$network), dsys$D, dsys$volume$cpp$xptr, dsys$T, trajectories)
sol2
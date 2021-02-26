TestPlotIndex <- function(){
    N <- 40
    DF <- PlotExample
    DFgrid <- expand.grid(X=0:N, Y=0:N)
    DFgrid <- dplyr::mutate(DFgrid, species=sample(c('A', 'B', 'C'), (N+1)^2, replace=TRUE),
        D_cm=runif((N+1)^2, 10, 50))
    DFgrid <- dplyr::mutate(DFgrid, H_m=10*sqrt(D_cm-5))
    Td <- TabDist(DFgrid, coord=c(0, N, 0, N))

    DFgrid2 <- dplyr::mutate(DFgrid, P=ceiling(3*(X+Y+1)/(2*(N+1))))
    DFgrid2 <- dplyr::mutate(DFgrid2, species=c('A','B','C')[P])
    Td2 <- TabDist(DFgrid2, coord=c(0, N, 0, N))

    DFgrid3 <- dplyr::mutate(DFgrid, species=rep(c('A', 'B'), ((N+2)^2)/2)[1:(N+1)^2])
    Td3 <- TabDist(DFgrid3, coord=c(0, N, 0, N))

    print(c('random', 'regular','irregular'))
    print(Compute_mingling(Td))
    print(Compute_mingling(Td2))
    print(Compute_mingling(Td3))

    DFgrid2 <- dplyr::mutate(DFgrid, P=ceiling(3*(X+Y+1)/(2*(N+1))))
    DFgrid2 <- dplyr::mutate(DFgrid2, D_cm=c(10, 20, 30)[P])
    Td2 <- TabDist(DFgrid2, coord=c(0, N, 0, N))

    DFgrid3 <- dplyr::mutate(DFgrid, D_cm=rep(seq(10,60,length.out=2),
	 ((N+2)^2)/2)[1:(N+1)^2])
    Td3 <- TabDist(DFgrid3, coord=c(0, N, 0, N))
    print(c('random', 'regular','irregular'))
    print(Compute_Size_Diff(Td))
    print(Compute_Size_Diff(Td2))
    print(Compute_Size_Diff(Td3))

    DFgrid2 <- dplyr::mutate(DFgrid, X=runif((N+1)^2, 0, N),
        Y=runif((N+1)^2, 0, N))
    Td2 <- TabDist(DFgrid2, coord=c(0, N, 0, N))

    nclust <- function(x0, y0, radius, n) {
        return(spatstat::runifdisc(n, radius, centre=c(x0, y0)))
    }
    PP <- spatstat::rPoissonCluster((N+5)^2/7, 0.2, nclust, radius=0.01, n=7)
    DFgrid3 <- dplyr::mutate(DFgrid, X=40*PP$x[1:41^2], Y=40*PP$y[1:41^2])
    Td3 <- TabDist(DFgrid3, coord=c(0, N, 0, N))
    print(c('grid', 'random','cluster'))
    print(Compute_Winkelmass(Td))
    print(Compute_Winkelmass(Td2))
    print(Compute_Winkelmass(Td3))
}

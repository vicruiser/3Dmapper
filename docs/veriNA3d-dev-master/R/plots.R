#' Plot a scatter&frequency circular plot for angular data
#'
#' For a vector of angular data (0 to 360), the function plots the 
#' distribution in a circular format.
#'
#' @param data A numeric vector with the data to plot.
#' @param clockwise A logical indicating the sense in which the data should be
#'     displayed.
#' @param start.degree An integer with the position in which the data starts
#'     being ploted.
#' @param main A string with the title of the plot.
#'
#' @return A circular plot with the input data.
#'
#' @examples
#'     ntinfo <- pipeNucData("1bau")
#'     C3endo_ntID <- cleanByPucker(ntinfo, pucker="C3'endo")
#'     C3endo <- ntinfo[ntinfo$ntID %in% C3endo_ntID, ]
#'     plotCircularDistribution(C3endo[, "delta"])
#'
#' @author Diego Gallego
#'
plotCircularDistribution <-
function(data, clockwise=FALSE, start.degree=0, main=NULL) {

    ## Prepare data to and labels to fit arguments ---------------------------
    if (!clockwise) {
        data <- abs(data - 360)
        labels <- append("", seq(from=330, to=0, by=-30))
    } else {
        labels <- append(seq(from=0, to=330, by=30), "")
    }

    ## Do the plot -----------------------------------------------------------
    fac <- factor(rep(1, times=length(data)))
    circos.par(start.degree=start.degree)
    circos.initialize(factors=fac, x=data, xlim=c(0, 360))
    circos.par("track.height"=0.05)
    circos.trackPlotRegion(factors=fac, ylim=c(0, 0.1), bg.border="white",
        panel.fun=function(x, y) {
            circos.axis(major.at=seq(from=0, to=360, by=30), labels=labels)
        })

    circos.trackPoints(fac, data, y=rep(0.05, times=length(data)),
                        col="blue", pch=16, cex=0.5)
    circos.par("track.height"=0.2)
    circos.trackHist(fac, data, bg.col="white", bg.border="white",
                        col=rgb(0.1, 0.5, 0.8, 0.3), breaks=360)

    ## Print title to plot ---------------------------------------------------
    if (!is.null(main)) {
        title(main=main)
    }

    circos.clear()
}
##############################################################################

#' Barplot wrapper
#'
#' Function to make more straigtforward the process of ploting a barplot for
#' categorical data.
#'
#' @param ntinfo A data.frame with the input data. It should contain the 
#'     columns with the desired categorical data and a column labeled ntID.
#' @param ntID A vector of integers with the desired nucleotides of 
#'     analysis. If NULL all the nucleotides in the data.frame will be used.
#' @param field The column name with the desired data.
#' @param na.rm A logical to remove missing data.
#' @param main A string with the title of the plot.
#' @param cex To be passed to the par() function.
#' @param file A string with the name of the output file. If NULL, the
#'     plot will be printed to screen.
#' @param width The width of the plot (passed to the png() function).
#' @param height The height of the plot (passed to the png() function).
#' @param bg The background color of the plot (passed to the png() function).
#' @param units The unit to measure height and width (passed to the png() 
#'     function).
#' @param res Resolution (passed to the png() function).
#'
#' @return A barplot with the categorical data of interest, which can be
#'     directly saved  to a ".png" file.
#'
#' @examples
#'     ## To see all the types of trinucleotides in the dataset:
#'     ntinfo <- pipeNucData("1bau")
#'     plotCategorical(ntinfo=ntinfo, field="localenv")
#' 
#' @author Diego Gallego
#'
plotCategorical <-
function(ntinfo, field, ntID=NULL, na.rm=FALSE,
            main=NULL, cex=0.5,
            file=NULL, width=15, height=15,
            bg="white", units="cm", res=200) {

    ## Subset the data of interest -------------------------------------------
    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    data <- ntinfo[which(ntinfo$ntID %in% ntID), field]

    ## Replace missing data by - or remove it --------------------------------
    if (sum(is.na(data)) > 0) {
        if (na.rm) {
            data <- data[complete.cases(data)]
        } else {
            data[is.na(data)] <- "-"
        }
    }

    ## Make the plot ---------------------------------------------------------
    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg, units=units, res=res)
    }
    .plot_hist(data=data, main=main, cex=cex)

    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################

#' Make 2D data plot or 3D view of the density
#'
#' Function to make a plot highlighting the most populated regions.
#'
#' @param ntinfo A data.frame with the input data. It should contain the 
#'     columns with the desired data and a column labeled ntID.
#' @param ntID A vector of integers with the desired nucleotides of 
#'     analysis. If NULL all the nucleotides in the data.frame will be used.
#' @param dens The output of a kernel density estimation (e.g. kde2d function)
#'     over the data of interest. If it is NULL and drawcontour=TRUE, it will 
#'     be computed on the fly.
#' @param bandwidths In case dens=NULL and drawcontour=TRUE, it will be passed
#'     to kde2d().
#' @param x A string with the parameter to be placed in the x axis.
#' @param y A string with the parameter to be placed in the y axis.
#' @param drawcontour A logical to highlight the most populated regions of the
#'     plot.
#' @param sd_over_mean_contours A numeric vector with the standard deviations
#'     over the mean to plot the contours (in case drawcontour=TRUE).
#' @param etatheta A logical to return a eta-theta plot.
#' @param points A vector of integers. It should contain the ntID of the 
#'     nucleotides to print in a different color.
#' @param colpoints A string with the desired color if points is not NULL.
#' @param defaultview A string to set different default options. Chose
#'     between '2Dupview', 'leftview' and 'rightview'.
#' @param thetaplot Argument to be pased to `persp3D()` to set the view 
#'     perspective.
#' @param phiplot Argument to be pased to `persp3D()` to est the view
#'     perspective.
#' @param cleanerview A logical to remove the lower density region to get a 
#'     cleaner plot.
#' @param file A string with the name of the output file. If NULL, the
#'     plot will be printed to screen.
#' @param width The width of the plot (passed to the png() function).
#' @param height The height of the plot (passed to the png() function).
#' @param bg The background color of the plot (passed to the png() function).
#' @param units The unit to measure height and width (passed to the png() 
#'     function).
#' @param res Resolution (passed to the png() function).
#' @param xlab String for X axis label
#' @param ylab String for Y axis label
#' @param zlab String for Z axis label
#'
#' @return A plot in screen, which can be directly saved  to a ".png" file.
#'     * {plot2D} A scatter plot.
#'     * {plot3Ddens} A density map of the data in 3D.
#'
#' @examples
#'     ntinfo <- pipeNucData("1bau")
#'     C3endo_ntID <- cleanByPucker(ntinfo, pucker="C3'endo")
#'     plot2D(ntinfo=ntinfo, ntID=C3endo_ntID)
#' 
#' @author Diego Gallego
#'
#' @name plot2D
NULL
##############################################################################

#' @export
#' @rdname plot2D
plot2D <-
function(ntinfo, ntID=NULL, dens=NULL, bandwidths=NULL, 
            x="eta", y="theta", drawcontour=TRUE, 
            sd_over_mean_contours=c(1, 2, 4), etatheta=FALSE,
            points=NULL, colpoints="red", 
            file=NULL, width=15, height=15,
            bg="white", units="cm", res=200, xlab="default", ylab="default") {

    ## Check input -----------------------------------------------------------
    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    if (!x %in% colnames(ntinfo) | !y %in% colnames(ntinfo)) {
        stop("Provide strings in the x&y arguments that match ", 
                "two columns of the data.frame ntinfo", sep="")
    }

    ## Save x and y values ---------------------------------------------------
    X <- x
    Y <- y
    inds <- which(complete.cases(ntinfo[, c(x, y)]))
    useful <- ntinfo[inds, "ntID"]
    ntID <- ntID[ntID %in% useful]
    x <- ntinfo[ntinfo$ntID %in% ntID, x]
    y <- ntinfo[ntinfo$ntID %in% ntID, y]

    ## If neccessary, find high-density regions ------------------------------
    if (drawcontour) {
        if (is.null(dens)) {
            if (is.null(bandwidths)) {
                bandwidths <- c(40, 40)
            }
            dens <- kde2d(x, y, 
                            n=c(361, 361), h=bandwidths, 
                            lims=c(0, 360, 0, 360))
        }
        mean_z <- mean(dens$z)
        sd_z=sd(dens$z)
    }

    ## Make the plot ---------------------------------------------------------
    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg, units=units, res=res)
    }

    if (etatheta) {
        if (xlab == "default") {
            xlab <- expression(paste(eta, " (degrees)", sep=""))
        }
        if (ylab == "default") {
            ylab <- expression(paste(theta, " (degrees)", sep=""))
        }
        plot(x, y,
                xlim=c(0, 360), ylim=c(0, 360), xlab=xlab, ylab=ylab,
                pch=19, cex=0.3, col="gray70", xaxt="n", yaxt="n")
    } else {
        if (xlab == "default") {
            xlab <- paste(X, " (degrees)", sep="")
        }
        if (ylab == "default") {
            ylab <- paste(Y, " (degrees)", sep="")
        }
        plot(x, y,
                xlim=c(0, 360), ylim=c(0, 360), xlab=xlab, ylab=ylab,
                pch=19, cex=0.3, col="gray70", xaxt="n", yaxt="n")
    }

    ## Add vertical and horizontal lines to highlight helical region ---------
    if (etatheta) {
        abline(h=190, lty=2, lwd=1.5, col="red", cex=2)
        abline(h=240, lty=2, lwd=1.5, col="red", cex=2)
        abline(v=150, lty=2, lwd=1.5, col="red", cex=2)
        abline(v=190, lty=2, lwd=1.5, col="red", cex=2)
    }

    ## Finish axis -----------------------------------------------------------
    axis(1, labels=seq(0, 360, by=36), at=seq(0, 360, by=36), las=2)
    axis(2, labels=seq(0, 360, by=36), at=seq(0, 360, by=36), las=2)

    ## If necessary, add points in a different color -------------------------
    if (!is.null(points)) {
        points <- points[points %in% useful]
        points(ntinfo[ntinfo$ntID %in% points, X],
                ntinfo[ntinfo$ntID %in% points, Y],
                col=colpoints, pch=19, cex=0.3)
    }

    ## Add contour lines to highlight high-density regions -------------------
    if (drawcontour) {
        levels <- mean_z + sd_over_mean_contours * sd_z
        if (length(sd_over_mean_contours) == 3 &&
                    sum(sd_over_mean_contours == c(1, 2, 4)) == 3) {

            colors <- c("cadetblue1", "deepskyblue", "blue4")
            legendtxt <- c(expression(paste("<", rho, ">+1s.d.")),
                            expression(paste("<", rho, ">+2s.d.")),
                            expression(paste("<", rho, ">+4s.d.")))
        } else {
            colors <- suppressWarnings(brewer.pal(length(levels), "Blues"))
            legendtxt <- paste("<mean>+", sd_over_mean_contours, "sd", sep="")
        }

        contour(dens, levels=levels, col=colors,
                add=TRUE, lwd=2, drawlabels=FALSE)
        legend("bottomleft", legend=legendtxt, 
                    lty=1, lwd=1, bty="n", col=colors)
    }

    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################

#' @export
#' @rdname plot2D
plot3Ddens <-
function(ntinfo, ntID=NULL, dens=NULL, bandwidths=NULL,
            x="eta", y="theta",
            defaultview=NULL, 
            thetaplot, phiplot, cleanerview=FALSE,
            file=NULL, width=15, height=15,
            bg="white", units="cm", res=600, 
            xlab="default", ylab="default", zlab="default") {

    ## Check input -----------------------------------------------------------
    ntID <- .giveMeValidntIDs(ntinfo, ntID)
    if (!x %in% colnames(ntinfo) | !y %in% colnames(ntinfo)) {
        stop("Provide strings in the x&y arguments that match ", 
                "two columns of the data.frame ntinfo", sep="")
    }

    ## Save x and y values ---------------------------------------------------
    inds <- which(complete.cases(ntinfo[, c(x, y)]))
    useful <- ntinfo[inds, "ntID"]
    ntID <- ntID[ntID %in% useful]
    x0 <- ntinfo[ntinfo$ntID %in% ntID, x]
    y0 <- ntinfo[ntinfo$ntID %in% ntID, y]

    ## If neccessary, find high-density regions ------------------------------
    if (is.null(dens)) {
        if (is.null(bandwidths)) {
            bandwidths <- c(40, 40)
        }
        dens <- kde2d(x0, y0, 
                        n=c(361, 361), h=bandwidths, 
                        lims=c(0, 360, 0, 360))
    }

    mean_z <- mean(dens$z)
    sd_z=sd(dens$z)

    ## Find range of data ----------------------------------------------------
    etarange <- range(x0)
    etaseq <- seq(etarange[1], etarange[2], length=361)
    thetarange <- range(y0)
    thetaseq <- seq(thetarange[1], thetarange[2], length=361)
    newdensZ <- dens$z

    ## Set parameters for defaultviews ---------------------------------------
    if (!is.null(defaultview)) {
        if (!defaultview %in% c("2Dupview", "leftview", "rightview"))
            stop("defaultview should be '2Dupview', 'leftview', 'rightview'")

        if (defaultview == "2Dupview") {
            cleanerview <- FALSE
            thetaplot <- 0
            phiplot <- 90
        } else if (defaultview == "leftview") {
            cleanerview <- TRUE
            thetaplot <- -30
            phiplot <- 30
        } else if (defaultview == "rightview") {
            cleanerview <- TRUE
            thetaplot <- 100
            phiplot <- 40
        }
    } else {
        if (missing(thetaplot))
            thetaplot <- -30
        if (missing(phiplot))
            phiplot <- 30
    }

    ## Remove low density data to have a cleaner view ------------------------
    if (cleanerview) {
        for (j in seq_len(ncol(newdensZ))) {
            for (i in seq_len(nrow(newdensZ))) {
                if (newdensZ[i, j] < mean_z) {
                    newdensZ[i, j] <- NA
                }
            }
        }
    }

    ## Make plot -------------------------------------------------------------
    if (!is.null(file)) {
        png(file, width=width, height=height, bg=bg, units=units, res=res)
    }

    if (xlab == "default") {
        xlab <- x
    }
    if (ylab == "default") {
        ylab <- y
    }
    if (zlab == "default") {
        zlab <- ""
    }
    par(mfrow=c(1, 1), mar=c(3, 3, 1.0, 3), cex=0.7, lty=1, las=1)
    persp3D(x=etaseq,
            y=thetaseq,
            z=newdensZ,
            border=NA, theta=thetaplot, phi=phiplot,
            xlab=xlab, ylab=ylab, zlab=zlab, lighting=TRUE)

    if (!is.null(file)) {
        dev.off()
    }
}
##############################################################################

#' Find 2D High Density Regions and return list of nucleotides
#'
#' Given a data.frame ("ntinfo") the function computes a 2D Kernel Density 
#' Estimation between the desired fields and returns a list of 
#' nucleotides (according with the ntID column) that 
#' are found in High Density Regions of the 2D map.
#'
#' @param ntID an object of class integer with the desired nucleotides of 
#'     analysis. If NULL all the nucleotides in the data.frame will be used.
#' @param ntinfo a data.frame with the input data. It should contain three 
#'     columns (additional columns will be ignored). One of them should be 
#'     "ntID" and the other two are optional and can be specified using the
#'     parameters "x" and "y".
#' @param x name of the column that will be used as "x" axis.
#' @param y name of the column that will be used as "y" axis.
#' @param SD_DENS height above the mean to be used to select the nucleotides
#' @param bandwidths object to be passed to the "kde2d" function (only used if
#'     "dens" is NULL).
#' @param dens optional object containing the output of "kde2d" or equivalent.
#' @param lims The limits of the rectangle covered by the grid as c(xl, xu,
#'     yl, yu).
#'
#' @return a list of vectors containing the nucleotides ID (according with the 
#'     ntID column) that clusterize in different regions of the 2D diagram.
#'
#' @author Diego Gallego
#'

#' @export
#' @rdname findHDR
findHDR <-
function(ntID=NULL, ntinfo, x="eta", y="theta", SD_DENS=1, 
            bandwidths=c(40, 40), dens=NULL, lims=c(0, 360, 0, 360)) {

    ## Make sure to use only valid nucleotides
    ntID <- .giveMeValidntIDs(ntinfo, ntID)

    ## Compute kernel density estimation if not provided
    if (is.null(dens)) {
        dens=kde2d(ntinfo[ntID, x], ntinfo[ntID, y],
                    n=c(length(lims[1]:lims[2]), length(lims[3]:lims[4])),
                    h=bandwidths, lims=lims)
    }
    mean_dens=mean(dens$z)
    sd_dens=sd(dens$z)

    ## Find the cells of the matrix "dens" with density above desired
    grid_cells <- which(dens$z>mean_dens + SD_DENS * sd_dens, arr.ind=TRUE)
    grid_cells <- grid_cells[order(grid_cells[, 2]),]
    grid_cells <- grid_cells[order(grid_cells[, 1]),]

    ## For the dense cells, find their dense neigbours
    dis_map <- nn2(grid_cells, k=9)
    cell_contacts <- c()
    for (i in seq_len(nrow(grid_cells))) {
        cell_contacts[[i]] <- dis_map$nn.idx[i, dis_map$nn.dists[i, ] < 2]
    }

    ## Create first cluster
    cluster <- 1
    clusters <- list()
#    clusters[[cluster]] <- paste(grid_cells[1, ], collapse="_")
    clusters[[cluster]] <- 1
    for (i in seq(2, nrow(grid_cells))) {
        point_contacts <- cell_contacts[[i]]

        ## Check if point pertains to an existing cluster
        for (j in seq_len(cluster)) {
            cluster_points <- unique(unlist(cell_contacts[clusters[[j]]]))
            if (any(point_contacts %in% cluster_points)) {
                clusters[[j]] <- append(clusters[[j]], i)
                assigned <- TRUE
                break()
            } else {
                assigned <- FALSE
            }
        }

        ## When the point is not directly connected to the existing clusters
        if (!assigned) {
            cluster <- cluster + 1
            clusters[[cluster]] <- i
        }
    }

    ## At this point, one same cluster might have been classified as two
    definitive <- vector(mode="integer", length=cluster)
    definitive[1] <- 1
    clusters2 <- list()
    clusters2[[1]] <- clusters[[1]]
    cluster2 <- 1

    for (i in seq_len(cluster)) {
        cluster_points <- unique(unlist(cell_contacts[clusters[[i]]]))

        ## If not assigned already, check if it is from a saved new cluster
        if (definitive[i] == 0) {
            ## Check if it pertains to previous saved clusters
            newclusts <- unique(definitive[definitive != 0])
            for (j in newclusts) {
                cluster_points2 <- unique(unlist(cell_contacts[clusters2[[j]]]))
            
                if (any(cluster_points %in% cluster_points2)) {
                    current_cluster <- definitive[j]
                    clusters2[[current_cluster]] <- 
                        append(clusters2[[current_cluster]], clusters[[i]])

                    definitive[i] <- current_cluster
                    break()
                }
            }

            ## If not assigned after this process, create new cluster
            if (definitive[i] == 0) {
                cluster2 <- cluster2 + 1
                clusters2[[cluster2]] <- clusters[[i]]
                definitive[i] <- cluster2
            }
        }

        ## In either case, see what else is in same cluster
        current_cluster <- definitive[i]

        for (j in seq_len(cluster)[-i]) {
            if (definitive[j] == 0) {
                cluster_points2 <- unique(unlist(cell_contacts[clusters[[j]]]))
                if (any(cluster_points %in% cluster_points2)) {
                    clusters2[[current_cluster]] <- 
                        append(clusters2[[current_cluster]], clusters[[j]])
                    definitive[j] <- i
                }
            }
        }
    }

    ## Assign definite variables
    cluster <- cluster2
    clusters <- clusters2

    ## Replace 1D indices for 2D indices
    for (i in seq_len(cluster)) {
        clusters[[i]] <- unlist(lapply(clusters[[i]], 
                                function(j) {
                                    paste(grid_cells[j, ], collapse="_")
                                }))
    }

    ## Find 2D indices for input angles
    angle_x_intervals <- seq(lims[1], lims[2], 
                                by=lims[2]/length(lims[1]:lims[2]))
    angle_y_intervals <- seq(lims[3], lims[4], 
                                by=lims[4]/length(lims[3]:lims[4]))

    ## Each of the points is assigned to one of the cells of the matrix
    grid <- unlist(lapply(ntID, 
                            FUN=function(.ntID) {
                                gridx <- which(ntinfo[ntinfo$ntID == .ntID, x] 
                                                < angle_x_intervals)[1] - 1
                                gridy <- which(ntinfo[ntinfo$ntID == .ntID, y] 
                                                < angle_y_intervals)[1] - 1
                                return(paste(gridx, gridy, sep="_"))
                            }))
    grid_coords <- as.data.frame(cbind(ntID, grid), stringsAsFactors=FALSE)
    grid_coords$ntID <- as.numeric(grid_coords$ntID)

    ## The cells found for each cluster are compared with the cells of each
    ## point and the population for each cluster is found
    output <- lapply(clusters, 
                        FUN=function(.cluster) {
                            ind <- which(grid_coords$grid %in% .cluster)
                            return(grid_coords[ind, "ntID"])
                        })
    return(output)
}
##############################################################################
## ===========================================================================
## Subfunctions
## ===========================================================================

.plot_hist <-
function(data, main=NULL, cex=0.5) {

    par(mfrow=c(1, 1))
    dataFactor <- as.factor(data)
    labels <- as.numeric(round(100 *
                                table(dataFactor)/sum(table(dataFactor)), 1))
    ylim=c(0, 1.1 * max(table(dataFactor)))

    xx <- barplot(table(dataFactor), main=main, 
                    density=TRUE, ylim=ylim, xaxt="n")

    text(xx, y=table(dataFactor), 
            labels=paste(labels, "%", sep=""), pos=3, cex=cex)
    axis(1, at=xx, labels=names(table(dataFactor)), 
            tick=FALSE, las=2, cex.axis=cex)
}
## ===========================================================================

.giveMeValidntIDs <-
function(ntinfo, ntID) {

    if (is.null(ntID)) {
        ntID <- ntinfo[, "ntID"]
    } else {
        if (sum(ntID %in% ntinfo$ntID)!=length(ntID)) {
            stop("Some of the specified IDs does not match an existing ID",
                    " in your data.frame", sep="")
        }
    }
    return(ntID)
}
## ===========================================================================

.onAttach <- function(lib, pkg) {
    msg <- paste("\nThanks for using veriNA3d!\n",
                "If you find it useful for your research, please cite us.\n",
                "Find the details typing \"citation('veriNA3d')\".\n", sep="")
    packageStartupMessage(msg)
}

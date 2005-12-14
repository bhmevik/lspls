.First.lib <- function(libname, pkgname) {
    ## FIXME: Create a separate options facility, like in 'sm'.
    if (is.null(getOption("pls.algorithm")))
        options(pls.algorithm = "oscorespls")
}

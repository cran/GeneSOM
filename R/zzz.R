.First.lib <- function(lib, pkg) {
  library.dynam( "GeneSOM", pkg, lib )
  require(mva)
}

.onLoad <- function(lib, pkg)
{
  library.dynam("passe", pkg, lib)
}

.onUnload <- function(lib)
{
  library.dynam.unload("passe", lib)
}

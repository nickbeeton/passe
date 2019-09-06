passe = function(pars, tree, n_species, n1 = 1L, n2 = 1L,
                 diag = FALSE, res.val = -Inf, fortran.func = 'passe')
{
  pars = abs(pars)
  n_pars = nrow(pars)
  if (is.null(n_pars))
  {
    n_pars = 1L
    pars = matrix(pars, 1, 6)
  }

  edge = tree$edge[ii <- order(tree$edge[,1]),]
  edge = as.integer(edge[,2])
  parent = integer(2*n_species - 1)
  parent[edge] = 1L:(2*n_species - 2L)
  parent = parent[-n_species:-1L]
  edge_length = tree$edge.length[ii]

  res = .Fortran(fortran.func,
                 l = numeric(n_pars),
                 edge = edge,
                 parent = parent,
                 edge_length = edge_length,
                 n_species = n_species,
                 traits = as.numeric(tree$tip.state),
                 n_splits1 = n1,
                 n_splits2 = n2,
                 mu0 = pars[,3],
                 mu1 = pars[,4],
                 q01 = pars[,5],
                 q10 = pars[,6],
                 l0 = pars[,1],
                 l1 = pars[,2],
                 n_pars = n_pars,
                 diag = diag
  )$l
  res[is.na(res)] = res.val
  res[is.infinite(res)] = res.val
  res
}

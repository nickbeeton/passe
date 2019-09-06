fast.bisse = function(n_species, pars, fortran_func = 'fast_bisse', diag = FALSE)
{
  l0 = pars[1]
  l1 = pars[2]
  mu0 = pars[3]
  mu1 = pars[4]
  q01 = pars[5]
  q10 = pars[6]

  g = l0 + mu1 - l1 - mu0
  if (g == 0.)
  {
    x = q10 / (q01 + q10)
  }  else  {
    b = g - q01 - q10
    x = (b + sqrt(b^2 + 4.*g*q10))/2./g
  }

  N.mean = exp(x*(l0-mu0) + (1-x)*(l1-mu1))

  const = (x*(l0+mu0+q01) + (1.-x)*(l1+mu1+q10))/(x*(l0-mu0) + (1.-x)*(l1-mu1))

  edges = 2L*n_species - 2L

  test = .Fortran(fortran_func,
                  x0 = 0.5,
                  n_max =  as.integer(10*const*n_species),
                  n_species = n_species,
                  edge = matrix(0L, edges, 2),
                  edge_length = numeric(edges),
                  edge_state = integer(edges),
                  node_state = integer(n_species - 1L),
                  tip_state = integer(n_species),
                  mu0 = pars[3],
                  mu1 = pars[4],
                  q01 = pars[5],
                  q10 = pars[6],
                  l0 = pars[1],
                  l1 = pars[2],
                  diag = diag)

  res = list(edge = test$edge,
             Nnode = n_species - 1,
             tip.label = paste('t', 1:n_species, sep = ''),
             node.label = paste('n', n_species + 1:(n_species - 1), sep=''),
             edge.length = test$edge_length,
             edge.state = test$edge_state,
             tip.state = test$tip_state,
             node.state = test$node_state)

  class(res) = 'phylo'
  res
}

#!/usr/bin/python
import bempp.api

# bempp Parameters
#print(bempp.api.global_parameters.hmat.coarsening)
bempp.api.global_parameters.hmat.eps = 1e-8

bempp.api.global_parameters.hmat.max_block_size = 2048
bempp.api.global_parameters.hmat.min_block_size = 21
bempp.api.global_parameters.hmat.max_rank = 30

bempp.api.global_parameters.quadrature.double_singular = 5 ###

bempp.api.global_parameters.quadrature.near.max_rel_dist = 2.0
bempp.api.global_parameters.quadrature.near.single_order = 5
bempp.api.global_parameters.quadrature.near.double_order = 5

bempp.api.global_parameters.quadrature.medium.max_rel_dist = 4.0
bempp.api.global_parameters.quadrature.medium.single_order = 5
bempp.api.global_parameters.quadrature.medium.double_order = 5

#bempp.api.global_parameters.quadrature.far.max_rel_dist = ?
bempp.api.global_parameters.quadrature.far.single_order = 5
bempp.api.global_parameters.quadrature.far.double_order = 5

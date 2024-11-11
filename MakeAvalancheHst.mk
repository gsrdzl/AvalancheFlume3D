avalanche_hst.cuda : \
	avalanche_hst.o \
	calculate_particle_parameters.o \
	cells.o \
	cell_dim.o \
	clear_forces.o \
	error_check.o \
	finish_Verlet.o \
	force_dep.o \
	force_flume.o \
	force_funnel.o \
	force_hop.o \
	force_part.o \
	force_walls.o \
	generate_particles.o \
	place_particles.o \
	get_accel.o \
	init_Verlet.o \
	rannew64.o \
	rotate_vectors.o \
	set_initial_velocity.o
	nvcc -O2 \
	avalanche_hst.o \
	calculate_particle_parameters.o \
	cells.o \
	cell_dim.o \
	clear_forces.o \
	error_check.o \
	finish_Verlet.o \
	force_dep.o \
	force_flume.o \
	force_funnel.o \
	force_hop.o \
	force_part.o \
	force_walls.o \
	generate_particles.o \
	get_accel.o \
	init_Verlet.o \
	place_particles.o \
	rannew64.o \
	rotate_vectors.o \
	set_initial_velocity.o \
	-o avalanche_hst.cuda -lm

avalanche_hst.o : avalanche_hst.cu cabecera.hu
	nvcc -c -O2  avalanche_hst.cu

calculate_particle_parameters.o : calculate_particle_parameters.cu cabecera.hu
	nvcc -c -O2  calculate_particle_parameters.cu

cells.o : cells.cu cabecera.hu
	nvcc -c -O2  cells.cu

cell_dim.o : cell_dim.cu cabecera.hu
	nvcc -c -O2  cell_dim.cu

clear_forces.o : clear_forces.cu cabecera.hu
	nvcc -c -O2  clear_forces.cu

error_check.o : error_check.cu cabecera.hu
	nvcc -c -O2  error_check.cu

finish_Verlet.o : finish_Verlet.cu cabecera.hu
	nvcc -c -O2  finish_Verlet.cu

force_dep.o : force_dep.cu cabecera.hu
	nvcc -c -O2  force_dep.cu

force_flume.o : force_flume.cu cabecera.hu
	nvcc -c -O2  force_flume.cu

force_funnel.o : force_funnel.cu cabecera.hu
	nvcc -c -O2  force_funnel.cu

force_hop.o : force_hop.cu cabecera.hu
	nvcc -c -O2  force_hop.cu

force_part.o : force_part.cu cabecera.hu
	nvcc -c -O2  force_part.cu

force_walls.o : force_walls.cu cabecera.hu
	nvcc -c -O2  force_walls.cu

generate_particles.o : generate_particles.cu cabecera.hu
	nvcc -c -O2  generate_particles.cu

get_accel.o : get_accel.cu cabecera.hu
	nvcc -c -O2  get_accel.cu

init_Verlet.o : init_Verlet.cu cabecera.hu
	nvcc -c -O2  init_Verlet.cu

place_particles.o : place_particles.cu cabecera.hu
	nvcc -c -O2  place_particles.cu

rannew64.o : rannew64.cu cabecera.hu
	nvcc -c -O2  rannew64.cu

rotate_vectors.o : rotate_vectors.cu cabecera.hu
	nvcc -c -O2  rotate_vectors.cu

set_initial_velocity.o : set_initial_velocity.cu cabecera.hu
	nvcc -c -O2  set_initial_velocity.cu

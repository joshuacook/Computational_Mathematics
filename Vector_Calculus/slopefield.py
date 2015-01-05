def slopefield(fn, tmin, tmax, dt, ymin, ymax, dy):
	"""Generator for the slopefield ticks"""
	t = tmin + 0.5 * dt
	while t < tmax:
		y = ymin + 0.5 * dy
		while y < ymax:
			# tick 
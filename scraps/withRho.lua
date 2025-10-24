local function withRhoSphericalToAbsoluteCartesianPosition(position)
	local r, theta, phi = vec3.components(position)
	local rho = math.sqrt(r ^ 2 + wormhole.throatRadius ^ 2)
	return vec3(
		rho * math.sin(theta) * math.cos(phi),
		rho * math.sin(theta) * math.sin(phi),
		rho * math.cos(theta)
	)
end

local function withRhoAlternateSphericalToAbsoluteCartesianPosition(position)
	local r, theta, phi = vec3.components(position)
	local rho = math.sqrt(r ^ 2 + wormhole.throatRadius ^ 2)
	return vec3(
		rho * math.sin(theta) * math.cos(phi),
		rho * math.cos(theta),
		rho * math.sin(theta) * math.sin(phi)
	)
end

local function rTangentToRhoTangent(position, tangent)
	local r = position.x
	local rho = math.sqrt(r ^ 2 + wormhole.throatRadius ^ 2)

	local dR = tangent.x

	local changedR = r + dR
	local changedRho = math.sqrt(changedR ^ 2 + wormhole.throatRadius)
	local dRho = changedRho - rho

	local dTheta, dPhi = tangent.y, tangent.z
	return vec3(dRho, dTheta, dPhi)
end

local function withRhoSphericalToAbsoluteCartesianTangent(position, tangent)
	local r, theta, phi = vec3.components(position)
	local dR, dTheta, dPhi = vec3.components(tangent)

	local rho = math.sqrt(r ^ 2 + wormhole.throatRadius ^ 2)
	local withRhoTangent = rTangentToRhoTangent(position, tangent)
	local dRho = withRhoTangent.x

	return vec3(
		-- Total differentials of withRhoSphericalToAbsoluteCartesianPosition's output
		(
			math.sin(theta) * math.cos(phi) * dRho +
			rho * math.cos(theta) * math.cos(phi) * dTheta -
			rho * math.sin(theta) * math.sin(phi) * dPhi
		),
		(
			math.sin(theta) * math.sin(phi) * dRho +
			rho * math.cos(theta) * math.sin(phi) * dTheta +
			rho * math.sin(theta) * math.cos(phi) * dPhi
		),
		(
			math.cos(theta) * dRho -
			rho * math.sin(theta) * dTheta
		)
	)
end

local function withRhoAlternateSphericalToAbsoluteCartesianTangent(position, tangent)
	local r, theta, phi = vec3.components(position)
	local dR, dTheta, dPhi = vec3.components(tangent)

	local rho = math.sqrt(r ^ 2 + wormhole.throatRadius ^ 2)
	local withRhoTangent = rTangentToRhoTangent(position, tangent)
	local dRho = withRhoTangent.x

	return vec3(
		-- Total differentials of withRhoAlternateSphericalToAbsoluteCartesianPosition's output
		(
			math.sin(theta) * math.cos(phi) * dRho +
			rho * math.cos(theta) * math.cos(phi) * dTheta -
			rho * math.sin(theta) * math.sin(phi) * dPhi
		),
		(
			math.cos(theta) * dRho -
			rho * math.sin(theta) * dTheta
		),
		(
			math.sin(theta) * math.sin(phi) * dRho +
			rho * math.cos(theta) * math.sin(phi) * dTheta +
			rho * math.sin(theta) * math.cos(phi) * dPhi
		)
	)
end

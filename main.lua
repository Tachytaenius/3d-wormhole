local mathsies = require("lib.mathsies")
local vec3 = mathsies.vec3
local quat = mathsies.quat

local consts = require("consts")

local wormhole
local camera
local sceneShader
local dummyTexture
local outputCanvas

local function euler(t, x, dt, f)
	x = x + f(t, x) * dt
	return x
end

local function rk4(t, x, dt, f)
	-- Not using divisions since the x value might not support it
	local k1 = f(t, x) * dt
	local k2 = f(t + dt * 0.5, x + k1 * 0.5) * dt
	local k3 = f(t + dt * 0.5, x + k2 * 0.5) * dt
	local k4 = f(t + dt, x + k3) * dt
	x = x + (k1 + k2 * 2 + k3 * 2 + k4) * (1 / 6)
	return x
end

local function sign(x)
	return x < 0 and -1 or x == 0 and 0 or 1
end

local function asinh(x)
	return math.log(x + math.sqrt(x ^ 2 + 1))
end

local function vectorLengthTangent(position, tangent)
	-- Line element of the metric tensor
	local r, theta, phi = vec3.components(position)
	local dR, dTheta, dPhi = vec3.components(tangent)
	return math.sqrt(dR ^ 2 + (wormhole.throatRadius ^ 2 + r ^ 2) * (dTheta ^ 2 + math.sin(theta) ^ 2 * dPhi ^ 2))
end

local function normaliseTangentVector(position, tangent)
	return tangent / vectorLengthTangent(position, tangent)
end

local function normaliseOrZeroTangent(position, tangent)
	local magnitude = vectorLengthTangent(position, tangent)
	if magnitude == 0 then
		return vec3() -- TODO: Is this right? It is in these cases but is it always? Tentatively, my intuition says that the answer has to do with positive-definiteness.
	end
	return tangent / magnitude
end

local function limitVectorLength(position, tangent, maxLength)
	local magnitude = vectorLengthTangent(position, tangent)
	if magnitude > maxLength then
		return normaliseTangentVector(position, tangent) * maxLength
	end
	return vec3.clone(tangent)
end

local function limitVectorLengthCartesian(v, maxLength)
	local magnitude = #v
	if magnitude > maxLength then
		return vec3.normalise(v) * maxLength
	end
	return vec3.clone(v)
end

-- I don't think these are strictly exactly spherical coordinates in this context but I am not sure what else to call them.
local function cartesianToSphericalBasis(position)
	local r, theta, phi = vec3.components(position)
	return
		vec3(1, 0, 0),
		vec3(0, 1 / math.sqrt(wormhole.throatRadius ^ 2 + r ^ 2), 0),
		vec3(0, 0, 1 / math.sqrt((wormhole.throatRadius ^ 2 + r ^ 2) * math.sin(theta) ^ 2))
end

local function sphericalToCartesianBasis(position)
	local r, theta, phi = vec3.components(position)
	return
		vec3(1, 0, 0),
		vec3(0, math.sqrt(wormhole.throatRadius ^ 2 + r ^ 2), 0),
		vec3(0, 0, math.sqrt((wormhole.throatRadius ^ 2 + r ^ 2) * math.sin(theta) ^ 2))
end

local function cartesianToSpherical(position, tangentCartesian)
	local xToR, yToTheta, zToPhi = cartesianToSphericalBasis(position)
	return tangentCartesian.x * xToR + tangentCartesian.y * yToTheta + tangentCartesian.z * zToPhi
end

local function sphericalToCartesian(position, tangentSpherical)
	local rToX, thetaToY, phiToZ = sphericalToCartesianBasis(position)
	return tangentSpherical.x * rToX + tangentSpherical.y * thetaToY + tangentSpherical.z * phiToZ
end

local function getChristoffelIndex(i, j, k)
	return i * 9 + j * 3 + k + 1
end

local function getChristoffelSymbols(position)
	-- Second kind

	local christoffels = {}
	for i = 1, 3 ^ 3 do
		christoffels[i] = 0
	end

	local r, theta, phi = vec3.components(position)

	christoffels[getChristoffelIndex(0, 1, 1)] = -r
	christoffels[getChristoffelIndex(0, 2, 2)] = -r * math.sin(theta) ^ 2

	christoffels[getChristoffelIndex(1, 1, 0)] = r / (wormhole.throatRadius ^ 2 + r ^ 2)
	christoffels[getChristoffelIndex(1, 0, 1)] = r / (wormhole.throatRadius ^ 2 + r ^ 2)
	christoffels[getChristoffelIndex(1, 2, 2)] = -math.sin(2 * theta) / 2

	christoffels[getChristoffelIndex(2, 0, 2)] = r / (wormhole.throatRadius ^ 2 + r ^ 2)
	christoffels[getChristoffelIndex(2, 2, 0)] = r / (wormhole.throatRadius ^ 2 + r ^ 2)
	christoffels[getChristoffelIndex(2, 1, 2)] = 1 / math.tan(theta)
	christoffels[getChristoffelIndex(2, 2, 1)] = 1 / math.tan(theta)

	return christoffels
end

local function getCameraRight()
	return cartesianToSpherical(camera.position,
		vec3.normalise(
			vec3.cross(
				sphericalToCartesian(camera.position, camera.forward),
				sphericalToCartesian(camera.position, camera.up)
			)
		)
	)
end

local stateMetatable
stateMetatable = {
	__add = function(a, b) -- b is expected to be another such array
		local new = setmetatable({}, stateMetatable)
		for i, av in ipairs(a) do
			new[i] = av + b[i]
		end
		return new
	end,
	__mul = function(a, b) -- b is expected to be a scalar
		local new = setmetatable({}, stateMetatable)
		for i, v in ipairs(a) do
			new[i] = v * b
		end
		return new
	end
}

function love.load()
	sceneShader = love.graphics.newShader("shaders/scene.glsl")
	dummyTexture = love.graphics.newImage(love.image.newImageData(1, 1))
	outputCanvas = love.graphics.newCanvas()

	wormhole = {
		throatRadius = 15
	}

	local position = vec3(-100, consts.tau / 4, 0)
	local forward = normaliseTangentVector(position, vec3(1, 0, 0))
	local up = normaliseTangentVector(position, vec3(0, 1, 0))
	camera = {
		position = position,
		velocity = vec3(),
		acceleration = 50,

		forward = forward,
		up = up,
		angularVelocity = vec3(),
		angularAcceleration = 1,

		verticalFOV = math.rad(90)
	}
end

function love.update(dt)
	local translation = vec3()
	local forward, up, right = camera.forward, camera.up, getCameraRight()
	if love.keyboard.isDown("d") then
		translation = translation + right
	end
	if love.keyboard.isDown("a") then
		translation = translation - right
	end
	if love.keyboard.isDown("e") then
		translation = translation + up
	end
	if love.keyboard.isDown("q") then
		translation = translation - up
	end
	if love.keyboard.isDown("w") then
		translation = translation + forward
	end
	if love.keyboard.isDown("s") then
		translation = translation - forward
	end
	local acceleration = normaliseOrZeroTangent(camera.position, translation) * camera.acceleration

	local rotation = vec3()
	if love.keyboard.isDown("k") then
		rotation = rotation + right
	end
	if love.keyboard.isDown("i") then
		rotation = rotation - right
	end
	if love.keyboard.isDown("l") then
		rotation = rotation + up
	end
	if love.keyboard.isDown("j") then
		rotation = rotation - up
	end
	if love.keyboard.isDown("u") then
		rotation = rotation + forward
	end
	if love.keyboard.isDown("o") then
		rotation = rotation - forward
	end
	rotation = -rotation -- TODO: Understand why this is needed when it isn't in other projects
	local angularAcceleration = limitVectorLength(camera.position, rotation, camera.angularAcceleration)

	local forwardCartesian = sphericalToCartesian(camera.position, forward)
	local upCartesian = sphericalToCartesian(camera.position, up)
	local rotationQuat = quat.fromAxisAngle(sphericalToCartesian(camera.position, angularAcceleration) * dt)
	-- local rotationQuat = quat.fromAxisAngle(angularAcceleration * dt)
	local forwardCartesianRotated = vec3.rotate(forwardCartesian, rotationQuat)
	local upCartesianRotated = vec3.rotate(upCartesian, rotationQuat)
	camera.forward = cartesianToSpherical(camera.position, forwardCartesianRotated)
	camera.up = cartesianToSpherical(camera.position, upCartesianRotated)

	camera.velocity = camera.velocity + acceleration * dt

	-- Move position and parallel transport tangent vectors
	local method = euler
	local steps = 1
	local stepSize = dt / steps
	-- Parallel transported tangent vectors (such as the forward vector) are updated separately after position and velocity as this seems to maintain accuracy
	local state = setmetatable({
		camera.position.x, camera.position.y, camera.position.z,
		camera.velocity.x, camera.velocity.y, camera.velocity.z
	}, stateMetatable)
	local parallelTransportState = setmetatable({
		camera.forward.x, camera.forward.y, camera.forward.z,
		camera.up.x, camera.up.y, camera.up.z
	}, stateMetatable)
	for i = 1, steps do
		local t = i * stepSize
		state = method(t, state, stepSize, function(t, state)
			local pos = vec3(state[1], state[2], state[3])
			local vel = vec3(state[4], state[5], state[6])
			local christoffels = getChristoffelSymbols(pos)

			local newState = setmetatable({}, stateMetatable)

			-- Position derivative
			newState[1] = vel.x
			newState[2] = vel.y
			newState[3] = vel.z

			-- Velocity derivative
			local accel = vec3()
			for i = 0, 2 do
				for j = 0, 2 do
					for k = 0, 2 do
						-- vec[consts.indexChars[index]] is just a way to number-index a vector that can normally only be indexed by letter. Note that x is 1, not 0.
						accel[consts.indexChars[i + 1]] = accel[consts.indexChars[i + 1]] -
							christoffels[getChristoffelIndex(i, j, k)] *
							vel[consts.indexChars[j + 1]] *
							vel[consts.indexChars[k + 1]]
					end
				end
			end
			newState[4] = accel.x
			newState[5] = accel.y
			newState[6] = accel.z

			return newState
		end)
		parallelTransportState = method(t, parallelTransportState, stepSize, function(t, parallelTransportState)
			local pos = vec3(state[1], state[2], state[3])
			local vel = vec3(state[4], state[5], state[6])
			local christoffels = getChristoffelSymbols(pos)
			local forward = vec3(parallelTransportState[1], parallelTransportState[2], parallelTransportState[3])
			local up = vec3(parallelTransportState[4], parallelTransportState[5], parallelTransportState[6])

			local newState = setmetatable({}, stateMetatable)

			local forwardDerivative = vec3()
			for i = 0, 2 do
				for j = 0, 2 do
					for k = 0, 2 do
						-- vec[consts.indexChars[index]] is just a way to number-index a vector that can normally only be indexed by letter. Note that x is 1, not 0.
						forwardDerivative[consts.indexChars[i + 1]] = forwardDerivative[consts.indexChars[i + 1]] -
							christoffels[getChristoffelIndex(i, j, k)] *
							forward[consts.indexChars[j + 1]] *
							vel[consts.indexChars[k + 1]]
					end
				end
			end
			newState[1], newState[2], newState[3] = vec3.components(forwardDerivative)

			local upDerivative = vec3()
			for i = 0, 2 do
				for j = 0, 2 do
					for k = 0, 2 do
						-- vec[consts.indexChars[index]] is just a way to number-index a vector that can normally only be indexed by letter. Note that x is 1, not 0.
						upDerivative[consts.indexChars[i + 1]] = upDerivative[consts.indexChars[i + 1]] -
							christoffels[getChristoffelIndex(i, j, k)] *
							up[consts.indexChars[j + 1]] *
							vel[consts.indexChars[k + 1]]
					end
				end
			end
			newState[4], newState[5], newState[6] = vec3.components(upDerivative)

			return newState
		end)
	end
	camera.position = vec3(state[1], state[2] % (consts.tau / 2), state[3] % consts.tau)
	camera.velocity = vec3(state[4], state[5], state[6])
	-- TODO: Normalise forward and up and also realign
	camera.forward = vec3(parallelTransportState[1], parallelTransportState[2], parallelTransportState[3])
	camera.up = vec3(parallelTransportState[4], parallelTransportState[5], parallelTransportState[6])



	-- TEMP
	-- do return end
	local cameraRight = getCameraRight()
	camera.up = normaliseTangentVector(camera.position,
		cartesianToSpherical(camera.position,
			vec3.normalise(
				vec3.cross(
					sphericalToCartesian(camera.position, cameraRight),
					sphericalToCartesian(camera.position, camera.forward)
				)
			)
		)
	)
	camera.up = normaliseTangentVector(camera.position, camera.up)
	camera.forward = normaliseTangentVector(camera.position, camera.forward)

	camera.position.y = camera.position.y % (consts.tau / 2)
	camera.position.z = camera.position.z % consts.tau
end

function love.draw()
	love.graphics.setShader(sceneShader)
	love.graphics.setCanvas(outputCanvas)
	love.graphics.clear()
	sceneShader:send("initialCameraPosition", {vec3.components(camera.position)})
	sceneShader:send("initialCameraForward", {vec3.components(
		sphericalToCartesian(camera.position, camera.forward)
	)})
	sceneShader:send("initialCameraUp", {vec3.components(
		sphericalToCartesian(camera.position, camera.up)
	)})
	sceneShader:send("initialCameraRight", {vec3.components(
		sphericalToCartesian(camera.position, getCameraRight())
	)})
	local screenAspectRatio = outputCanvas:getWidth() / outputCanvas:getHeight()
	sceneShader:send("cameraHorizontalDirectionExtent", math.tan(camera.verticalFOV / 2))
	sceneShader:send("cameraVerticalDirectionExtent", math.tan(camera.verticalFOV / 2) / screenAspectRatio)
	sceneShader:send("wormholeThroatRadius", wormhole.throatRadius)
	sceneShader:send("raySpeed", consts.rayStepSize)
	sceneShader:send("rayStepCount", consts.rayStepCount)
	love.graphics.draw(dummyTexture, 0, 0, 0, outputCanvas:getDimensions())
	love.graphics.setShader()
	love.graphics.setCanvas()

	love.graphics.draw(outputCanvas)

	love.graphics.print(love.timer.getFPS())
end

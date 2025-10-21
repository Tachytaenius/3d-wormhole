local mathsies = require("lib.mathsies")
local vec3 = mathsies.vec3
local quat = mathsies.quat
local mat4 = mathsies.mat4

local consts = require("consts")

local wormhole
local camera
local sceneShader
local dummyTexture
local outputCanvas
local upperSkybox, lowerSkybox

local renderModeInfo

local rayStepSize
local rayStepCount

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

local function safeAcos(x) -- For when you're sure the input is in [-1, 1]
	return math.acos(math.max(-1, math.min(1, x)))
end

local function fixorder(t)
	-- local function index(i, j)
	-- 	return i + 3 * j + 1
	-- end
	-- for i = 0, 3 do
	-- 	for j = 0, 3 do
	-- 		if i > j then
	-- 			t[index(i, j)], t[index(j, i)] = t[index(j, i)], t[index(i, j)]
	-- 		end
	-- 	end
	-- end
	-- return unpack(t)
	return mat4.components(mat4.transpose(mat4(unpack(t))))
end

local function swapcols2and3(t)
	local function index(i, j)
		return i + 3 * j + 1
	end
	for i = 0, 3 do
		t[index(1, i)], t[index(2, i)] = t[index(2, i)], t[index(1, i)]
	end
	return unpack(t)
end

local function sphericalJacobian(position) -- Returns mat4 because there are no mat3s in mathsies
	local r, theta, phi = vec3.components(position)
	return mat4(
		fixorder{math.sin(theta) * math.cos(phi),
		math.sin(theta) * math.sin(phi),
		math.cos(theta),
		0,

		r * math.cos(theta) * math.cos(phi),
		r * math.cos(theta) * math.sin(phi),
		-r * math.sin(theta),
		0,

		-r * math.sin(theta) * math.sin(phi),
		r * math.sin(theta) * math.cos(phi),
		0,
		0,

		0,
		0,
		0,
		1}
	)
end

local function alternateSphericalInverseJacobian(position)
	local r, theta, phi = vec3.components(position)
	return mat4(
		fixorder{math.sin(theta) * math.cos(phi),
		math.cos(theta) * math.cos(phi) / r,
		-math.sin(phi) / (math.sin(theta) * r),
		0,

		math.cos(theta),
		-math.sin(theta) / r,
		0,
		0,

		math.sin(theta) * math.sin(phi),
		math.cos(theta) * math.sin(phi) / r,
		math.cos(phi) / (math.sin(theta) * r),
		0,

		0,
		0,
		0,
		1}
	)
end

local function coordSwitchPosition(position)
	local r, theta, phi = vec3.components(position)
	local angleDirection = vec3(
		math.sin(theta) * math.cos(phi),
		math.sin(theta) * math.sin(phi),
		math.cos(theta)
	)
	return vec3(
		r,
		safeAcos(angleDirection.y),
		sign(angleDirection.z) * safeAcos(angleDirection.x / math.sqrt(angleDirection.x ^ 2 + angleDirection.z ^ 2))
	)
end

local function coordSwitchTangent(position, tangent)
	return
		alternateSphericalInverseJacobian(coordSwitchPosition(position)) *
		sphericalJacobian(position) *
		tangent
end

local function sphericalToAbsoluteCartesianPosition(position)
	local r, theta, phi = vec3.components(position)
	return vec3(
		r * math.sin(theta) * math.cos(phi),
		r * math.sin(theta) * math.sin(phi),
		r * math.cos(theta)
	)
end

local function alternateSphericalToAbsoluteCartesianPosition(position)
	local r, theta, phi = vec3.components(position)
	return vec3(
		r * math.sin(theta) * math.cos(phi),
		r * math.cos(theta),
		r * math.sin(theta) * math.sin(phi)
	)
end

local function switchCameraCoords()
	-- We use a multiplier in some code in update and draw based on whether the camera is in alternative coords mode.
	-- This may be a bad way to do it, since maybe we just need to switch the coords of tangent vectors differently...?
	-- TODO: Investigate!

	-- local initialPosition = camera.position

	local initialR = camera.position.x
	local fakePosition = vec3(1, camera.position.y, camera.position.z)

	local switchedPosition = coordSwitchPosition(camera.position)

	-- local prevxyz = camera.altCoords and alternateSphericalToAbsoluteCartesianPosition(camera.position) or sphericalToAbsoluteCartesianPosition(camera.position)

	-- local prevr = getCameraRight()
	local function handleTangent(name, lengthMode)
		local initialLen = vectorLengthTangent(camera.position, camera[name])
		local initialDr = camera[name].x
		camera[name] = coordSwitchTangent(fakePosition, camera[name])
		camera[name].x = initialDr
		if
			lengthMode == "normalise" or
			lengthMode == "maintain" and vectorLengthTangent(switchedPosition, camera[name]) > 0
		then
			camera[name] = normaliseTangentVector(switchedPosition, camera[name])
		end
		if lengthMode == "maintain" then
			camera[name] = camera[name] * initialLen
		end
		-- local newLen = vectorLengthTangent(switchedPosition, camera[name])
		-- print(name, "lendiff", newLen - initialLen)
	end
	handleTangent("velocity")
	handleTangent("angularVelocity")
	handleTangent("forward", "normalise")
	handleTangent("up", "normalise")

	camera.position = switchedPosition
	camera.position.x = initialR

	camera.altCoords = not camera.altCoords

	-- local currr = getCameraRight()
	-- print("rightlendiff",
	-- 	vectorLengthTangent(camera.position, currr) -
	-- 	vectorLengthTangent(initialPosition, prevr)
	-- )

	-- local currxyz = camera.altCoords and alternateSphericalToAbsoluteCartesianPosition(camera.position) or sphericalToAbsoluteCartesianPosition(camera.position)

	-- print("posdist", vec3.distance(prevxyz, currxyz))
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

local function init(parameters)
	parameters = parameters or {}

	sceneShader = love.graphics.newShader("shaders/scene.glsl")
	dummyTexture = love.graphics.newImage(love.image.newImageData(1, 1))
	outputCanvas = love.graphics.newCanvas(parameters.canvasWidth, parameters.canvasHeight) -- If no canvas size is specified, this is the window size

	local lowerPath = "backgrounds/lower/"
	lowerSkybox = love.graphics.newCubeImage({
		lowerPath .. "1.png",
		lowerPath .. "2.png",
		lowerPath .. "3.png",
		lowerPath .. "4.png",
		lowerPath .. "5.png",
		lowerPath .. "6.png"
	})
	local upperPath = "backgrounds/upper/"
	upperSkybox = love.graphics.newCubeImage({
		upperPath .. "1.png",
		upperPath .. "2.png",
		upperPath .. "3.png",
		upperPath .. "4.png",
		upperPath .. "5.png",
		upperPath .. "6.png"
	})

	wormhole = {
		throatRadius = 15
	}

	local position = vec3(-75, consts.tau / 4, 0)
	local forward = normaliseTangentVector(position, vec3(1, 0, 0))
	local up = normaliseTangentVector(position, vec3(0, 1, 0))
	camera = {
		altCoords = false,

		position = position,
		velocity = parameters.initialVelocity or vec3(),
		acceleration = 50,
		maxSpeed = 50,

		forward = forward,
		up = up,
		angularVelocity = parameters.initialAngularVelocity or vec3(),
		angularAcceleration = 1,
		maxAngularSpeed = 1,

		verticalFOV = math.rad(90)
	}
end

local function update(dt, rendering)
	if
		camera.position.y < consts.altCoordsProportion * consts.tau / 2 or
		camera.position.y > (1 - consts.altCoordsProportion) * consts.tau / 2
	then
		switchCameraCoords()
		-- print("Switched")
	end

	local multiplier = camera.altCoords and -1 or 1
	local forward, up, right = camera.forward, camera.up, multiplier * getCameraRight()

	local translation = vec3()
	if not rendering then
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
	end
	local acceleration = normaliseOrZeroTangent(camera.position, translation) * camera.acceleration

	local rotation = vec3()
	if not rendering then
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
	end
	rotation = -rotation -- TODO: Understand why this is there when other projects don't need it.
	if not love.keyboard.isDown("lctrl") then
		-- Accelerate and limit angular velocity
		local angularAcceleration = limitVectorLength(camera.position, rotation, camera.angularAcceleration)
		camera.angularVelocity = camera.angularVelocity + angularAcceleration * dt
		camera.angularVelocity = limitVectorLength(camera.position, camera.angularVelocity, camera.maxAngularSpeed)
	else
		-- Brake
		local angularSpeed = vectorLengthTangent(camera.position, camera.angularVelocity)
		angularSpeed = math.max(0, angularSpeed - camera.angularAcceleration * dt)
		camera.angularVelocity = angularSpeed * normaliseOrZeroTangent(camera.position, camera.angularVelocity)
	end

	local forwardCartesian = sphericalToCartesian(camera.position, forward)
	local upCartesian = sphericalToCartesian(camera.position, up)
	local rotationQuat = quat.fromAxisAngle(sphericalToCartesian(camera.position, multiplier * camera.angularVelocity) * dt)
	local forwardCartesianRotated = vec3.rotate(forwardCartesian, rotationQuat)
	local upCartesianRotated = vec3.rotate(upCartesian, rotationQuat)
	camera.forward = cartesianToSpherical(camera.position, forwardCartesianRotated)
	camera.up = cartesianToSpherical(camera.position, upCartesianRotated)

	if not love.keyboard.isDown("lshift") then
		-- Accelerate and then limit velocity
		camera.velocity = camera.velocity + acceleration * dt
		camera.velocity = limitVectorLength(camera.position, camera.velocity, camera.maxSpeed)
	else
		-- Brake
		local speed = vectorLengthTangent(camera.position, camera.velocity)
		speed = math.max(0, speed - camera.acceleration * dt)
		camera.velocity = speed * normaliseOrZeroTangent(camera.position, camera.velocity)
	end

	-- Move position and parallel transport tangent vectors
	local method = rk4
	local steps = 1
	local stepSize = dt / steps
	-- Parallel transported tangent vectors (such as the forward vector) are updated separately after position and velocity as this seems to maintain accuracy
	local state = setmetatable({
		camera.position.x, camera.position.y, camera.position.z,
		camera.velocity.x, camera.velocity.y, camera.velocity.z
	}, stateMetatable)
	local parallelTransportState = setmetatable({
		camera.forward.x, camera.forward.y, camera.forward.z,
		camera.up.x, camera.up.y, camera.up.z,
		camera.angularVelocity.x, camera.angularVelocity.y, camera.angularVelocity.z
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
			local angularVelocity = vec3(parallelTransportState[7], parallelTransportState[8], parallelTransportState[9])

			local newState = setmetatable({}, stateMetatable)

			local function move(variable, startIndex)
				local derivative = vec3()
				for i = 0, 2 do
					for j = 0, 2 do
						for k = 0, 2 do
							-- vec[consts.indexChars[index]] is just a way to number-index a vector that can normally only be indexed by letter. Note that x is 1, not 0.
							derivative[consts.indexChars[i + 1]] = derivative[consts.indexChars[i + 1]] -
								christoffels[getChristoffelIndex(i, j, k)] *
								variable[consts.indexChars[j + 1]] *
								vel[consts.indexChars[k + 1]]
						end
					end
				end
				for i = 1, 3 do
					newState[startIndex + i - 1] = derivative[consts.indexChars[i]]
				end
			end

			move(forward, 1)
			move(up, 4)
			move(angularVelocity, 7)

			return newState
		end)
	end
	camera.position = vec3(state[1], state[2], state[3])
	camera.velocity = vec3(state[4], state[5], state[6])
	camera.forward = vec3(parallelTransportState[1], parallelTransportState[2], parallelTransportState[3])
	camera.up = vec3(parallelTransportState[4], parallelTransportState[5], parallelTransportState[6])
	camera.angularVelocity = vec3(parallelTransportState[7], parallelTransportState[8], parallelTransportState[9])

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
	-- camera.up = normaliseTangentVector(camera.position, camera.up)
	camera.forward = normaliseTangentVector(camera.position, camera.forward)

	-- theta (camera.position.y) won't go out of bounds (barring absurd speeds, I guess) because we switch coordinates before it gets close to 0 or tau / 2
	camera.position.z = camera.position.z % consts.tau
end

local function draw()
	local multiplier = camera.altCoords and -1 or 1

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
		sphericalToCartesian(camera.position, multiplier * getCameraRight())
	)})
	local screenAspectRatio = outputCanvas:getWidth() / outputCanvas:getHeight()
	sceneShader:send("cameraHorizontalDirectionExtent", math.tan(camera.verticalFOV / 2))
	sceneShader:send("cameraVerticalDirectionExtent", math.tan(camera.verticalFOV / 2) / screenAspectRatio)
	sceneShader:send("wormholeThroatRadius", wormhole.throatRadius)
	sceneShader:send("raySpeed", rayStepSize)
	sceneShader:send("rayStepCount", rayStepCount)
	sceneShader:send("initialAltCoords", camera.altCoords)
	sceneShader:send("altCoordsProportion", consts.altCoordsProportion)
	sceneShader:send("lowerSkybox", lowerSkybox)
	sceneShader:send("upperSkybox", upperSkybox)
	love.graphics.draw(dummyTexture, 0, 0, 0, outputCanvas:getDimensions())
	love.graphics.setShader()
	love.graphics.setCanvas()
end

function love.load(args)
	if args[1] == "--render" then
		rayStepSize = 0.25
		rayStepCount = 512

		local canvasWidth = 3840 -- 1920
		local canvasHeight = 2160 -- 1080
		local time = 20
		init({
			canvasWidth = canvasWidth,
			canvasHeight = canvasHeight,
			initialVelocity = vec3(80 / time, 0.1 / time, 0),
			initialAngularVelocity = vec3(0, 0.002 / time, 0)
		})
		local framerate = 60
		renderModeInfo = {
			ticksDone = 0,
			ticksRequired = math.ceil(time * framerate),
			tickDt = 1 / framerate,
			preview = true
		}

		-- It doesn't go into the love save directory
		-- local extension = "mkv" -- "mp4"
		-- local dateTime = os.date("%Y-%m-%d %H-%M-%S")
		-- local outputFilename
		-- local i = 1
		-- repeat
		-- 	outputFilename = dateTime .. " render " .. i .. "." .. extension
		-- until not love.filesystem.getInfo(outputFilename)
		-- local outputFilenameSpaceFix = outputFilename:gsub(" ", "\\ ")

		local outPath = args[2]

		local ffmpegFile, errorMessage = io.popen(
			"ffmpeg -f rawvideo -r " .. framerate .. " -pix_fmt rgba -s " .. canvasWidth .. "x" .. canvasHeight .. " -i - -c:v libx264 -crf:v 23 -preset:v veryslow -filter:v format=yuv420p " .. outPath,
			love.system.getOS() == "Windows" and "wb" or "w" -- Apparently windows needs "wb"
		)
		if not ffmpegFile then
			error(errorMessage)
		end
		renderModeInfo.ffmpegFile = ffmpegFile
	else
		rayStepSize = 1
		rayStepCount = 128

		init()
	end
end

function love.update(dt)
	if renderModeInfo then
		if not renderModeInfo.finished then
			update(renderModeInfo.tickDt, true)
			draw()
			local imageData = love.graphics.readbackTexture(outputCanvas)
			renderModeInfo.ffmpegFile:write(imageData:getString())
			imageData:release()
			renderModeInfo.ticksDone = renderModeInfo.ticksDone + 1
			if renderModeInfo.ticksDone >= renderModeInfo.ticksRequired then
				-- Finish
				renderModeInfo.ffmpegFile:close()
				renderModeInfo.finished = true
			end
		end
	else
		update(dt)
	end
end

function love.draw()
	if renderModeInfo then
		if renderModeInfo.preview then
			local cw, ch = outputCanvas:getDimensions() -- Canvas dimensions
			local sw, sh = love.graphics.getDimensions() -- Screen dimensions
			local s = math.min(sw / cw, sh / ch) -- Scale
			local x, y = (sw - cw * s) / 2, (sh - ch * s) / 2 -- Position
			love.graphics.draw(outputCanvas, x, y, 0, s, s)
		end
		love.graphics.print(
			"Tick progress: " .. renderModeInfo.ticksDone .. "/" .. renderModeInfo.ticksRequired .. "\n" ..
			(renderModeInfo.finished and "Finished!\n" or "")
		)
	else
		draw()
		love.graphics.draw(outputCanvas)
		love.graphics.print(love.timer.getFPS())
	end
end

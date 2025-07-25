local vec3 = require("lib.mathsies").vec3

local consts = {}

consts.tau = math.pi * 2

consts.rayStepSize = 0.5
consts.rayStepCount = 256

consts.forwardVector = vec3(0, 0, 1)
consts.upVector = vec3(0, 1, 0)
consts.rightVector = vec3(1, 0, 0)

consts.indexChars = {"x", "y", "z", "w"}

consts.poleAlternateChartAngle = consts.tau / 8  -- At this angular distance from the poles, switch to an alternate chart to avoid singularities

return consts

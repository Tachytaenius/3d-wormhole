local vec3 = require("lib.mathsies").vec3

local consts = {}

consts.tau = math.pi * 2

consts.forwardVector = vec3(0, 0, 1)
consts.upVector = vec3(0, 1, 0)
consts.rightVector = vec3(1, 0, 0)

consts.indexChars = {"x", "y", "z", "w"}

consts.poleAlternateChartAngle = consts.tau / 8  -- At this angular distance from the poles, switch to an alternate chart to avoid singularities
consts.altCoordsProportion = 1 / 16

-- Determines the distance(s) from the wormhole when the system switches between curved/flat mode
consts.wormholeCutoffGradient = 0.1
consts.wormholeCutoffExtraFactor = 0.1

return consts

const float tau = 6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413; // Thanks OEIS

uniform vec3 initialCameraForward;
uniform vec3 initialCameraUp;
uniform vec3 initialCameraRight;

#ifdef VERTEX

uniform float cameraHorizontalDirectionExtent;
uniform float cameraVerticalDirectionExtent;

layout(location = 0) in vec4 VertexPosition;

out vec3 directionPreNormalise;

void vertexmain() {
	vec4 transformedPosition = TransformProjectionMatrix * VertexPosition;
	// No w divide needed
	directionPreNormalise =
		initialCameraRight * cameraHorizontalDirectionExtent * transformedPosition.x +
		initialCameraUp * cameraVerticalDirectionExtent * transformedPosition.y +
		initialCameraForward;
	gl_Position = transformedPosition;
}

#endif

#ifdef PIXEL

in vec3 directionPreNormalise;

out vec4 outColour;

uniform int rayStepCount;
uniform float raySpeed;

uniform float wormholeThroatRadius;
uniform vec3 mouthAPosition;
uniform vec3 mouthBPosition;
uniform float wormholeCutoffGradient;
uniform float wormholeCutoffExtraFactor;

uniform vec3 initialCameraPosition;
uniform bool startInCurvedMode;
uniform bool startInUpperWorld;
uniform bool initialAltCoords;
uniform float altCoordsProportion;

uniform samplerCube lowerSkybox;
uniform samplerCube upperSkybox;

float safeAcos(float x) { // For use when, mathematically, the output should be within [-1, 1]
	return acos(clamp(x, -1.0, 1.0));
}

vec3 fixPosition(float r, vec3 converted) {
	if (r >= 0.0) {
		return mouthAPosition + converted;
	} else {
		vec3 delta = mouthBPosition - mouthAPosition;
		vec3 direction = normalize(delta);
		vec3 parallel = direction * dot(converted, direction);
		vec3 perpendicular = converted - parallel;
		vec3 parallelFlipped = -parallel;
		vec3 convertedFlipped = parallelFlipped + perpendicular;
		return mouthBPosition + convertedFlipped;
	}
}

// Reference frames:
// Relative cartesian (xyz aligned to change in r theta phi)
// Main spherical (r theta phi)
// Alternate spherical (r theta phi but the poles are perpendicular)
// Absolute cartesian (xyz as acquried from main spherical as though this were Euclidean 3D (so this one doesn't work unless far from the wormhole))

mat3 cartesianToSphericalBasis(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	float commonVar = wormholeThroatRadius * wormholeThroatRadius + r * r;
	float sine = sin(theta);
	return mat3(
		vec3(1.0, 0.0, 0.0),
		vec3(0.0, 1.0 / sqrt(commonVar), 0.0),
		vec3(0.0, 0.0, 1.0 / sqrt(commonVar * sine * sine))
	);
}

mat3 sphericalToCartesianBasis(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	float commonVar = wormholeThroatRadius * wormholeThroatRadius + r * r;
	float sine = sin(theta);
	return mat3(
		vec3(1.0, 0.0, 0.0),
		vec3(0.0, sqrt(commonVar), 0.0),
		vec3(0.0, 0.0, sqrt(commonVar * sine * sine))
	);
}

vec3 spherical1ToAbsoluteCartesianPosition(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	vec3 converted = vec3(
		r * sin(theta) * cos(phi),
		r * sin(theta) * sin(phi),
		r * cos(theta)
	);
	return fixPosition(r, converted);
}

vec3 spherical1ToAbsoluteCartesianTangent(vec3 position, vec3 tangent) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;

	float dR = tangent.x;
	float dTheta = tangent.y;
	float dPhi = tangent.z;

	return vec3(
		(
			sin(theta) * cos(phi) * dR +
			r * cos(theta) * cos(phi) * dTheta -
			r * sin(theta) * sin(phi) * dPhi
		),
		(
			sin(theta) * sin(phi) * dR +
			r * cos(theta) * sin(phi) * dTheta +
			r * sin(theta) * cos(phi) * dPhi
		),
		(
			cos(theta) * dR -
			r * sin(theta) * dTheta
		)
	);
}

vec3 spherical2ToAbsoluteCartesianPosition(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	vec3 converted = vec3(
		r * sin(theta) * cos(phi),
		r * cos(theta),
		r * sin(theta) * sin(phi)
	);
	return fixPosition(r, converted);
}

vec3 spherical2ToAbsoluteCartesianTangent(vec3 position, vec3 tangent) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;

	float dR = tangent.x;
	float dTheta = tangent.y;
	float dPhi = tangent.z;

	return vec3(
		(
			sin(theta) * cos(phi) * dR +
			r * cos(theta) * cos(phi) * dTheta -
			r * sin(theta) * sin(phi) * dPhi
		),
		(
			cos(theta) * dR -
			r * sin(theta) * dTheta
		),
		(
			sin(theta) * sin(phi) * dR +
			r * cos(theta) * sin(phi) * dTheta +
			r * sin(theta) * cos(phi) * dPhi
		)
	);
}

vec3 spherical1ToSpherical2Position(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	vec3 angleDirection = vec3(
		sin(theta) * cos(phi),
		sin(theta) * sin(phi),
		cos(theta)
	);
	return vec3(
		r,
		safeAcos(angleDirection.y),
		sign(angleDirection.z) * safeAcos(angleDirection.x / length(angleDirection.xz))
	);
}

vec3 spherical2ToSpherical1Position(vec3 position) {
	// The maths is the same
	return spherical1ToSpherical2Position(position);
}

mat3 spherical1Jacobian(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	return mat3(
		vec3(
			sin(theta) * cos(phi),
			sin(theta) * sin(phi),
			cos(theta)
		),
		vec3(
			r * cos(theta) * cos(phi),
			r * cos(theta) * sin(phi),
			-r * sin(theta)
		),
		vec3(
			-r * sin(theta) * sin(phi),
			r * sin(theta) * cos(phi),
			0.0
		)
	);
}

mat3 spherical1InverseJacobian(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;
	return mat3(
		vec3(
			sin(theta) * cos(phi),
			cos(theta) * cos(phi) / r,
			-sin(phi) / (sin(theta) * r)
		),
		vec3(
			sin(theta) * sin(phi),
			cos(theta) * sin(phi) / r,
			cos(phi) / (sin(theta) * r)
		),
		vec3(
			cos(theta),
			-sin(theta) / r,
			0.0
		)
	);
}

mat3 spherical2Jacobian(vec3 position) {
	mat3 for1 = spherical1Jacobian(position);
	return mat3(
		for1[0].xzy,
		for1[1].xzy,
		for1[2].xzy
	);
}

mat3 spherical2InverseJacobian(vec3 position) {
	mat3 for1 = spherical1InverseJacobian(position);
	return mat3(
		for1[0],
		for1[2],
		for1[1]
	);
}

vec3 spherical1ToSpherical2Velocity(vec3 position, vec3 velocity) {
	mat3 jacobian1 = spherical1Jacobian(position);
	mat3 inverseJacobian2 = spherical2InverseJacobian(spherical1ToSpherical2Position(position));
	return inverseJacobian2 * jacobian1 * velocity;
}

vec3 spherical2ToSpherical1Velocity(vec3 position, vec3 velocity) {
	return spherical1ToSpherical2Velocity(position, velocity);
}

// Coordinate system functions end

float getWormholeCutoffR(bool inside) {
	float cutoffR = (wormholeThroatRadius * sqrt(1.0 - wormholeCutoffGradient * wormholeCutoffGradient)) / wormholeCutoffGradient;

	if (inside) {
		return cutoffR * (1.0 + wormholeCutoffExtraFactor);
	} else {
		return cutoffR * (1.0 - wormholeCutoffExtraFactor);
	}
}

mat3[3] getChristoffelSymbols(vec3 position) {
	float r = position.x;
	float theta = position.y;
	float phi = position.z;

	float sinTheta = sin(theta);
	float cotTheta = 1.0 / tan(theta);
	float commonVar = r / (wormholeThroatRadius * wormholeThroatRadius + r * r);

	mat3 forR = mat3(0.0);
	forR[1][1] = -r;
	forR[2][2] = -r * sinTheta * sinTheta;

	mat3 forTheta = mat3(0.0);
	forTheta[1][0] = commonVar;
	forTheta[0][1] = commonVar;
	forTheta[2][2] = -sin(2.0 * theta) / 2.0;

	mat3 forPhi = mat3(0.0); 
	forPhi[0][2] = commonVar;
	forPhi[2][0] = commonVar;
	forPhi[1][2] = cotTheta;
	forPhi[2][1] = cotTheta;

	return mat3[3] (forR, forTheta, forPhi);
}

mat2x3 integrationStateDeriv(float t, mat2x3 state) {
	vec3 position = state[0];
	vec3 velocity = state[1];

	vec3 positionDerivative = velocity;

	mat3[3] christoffelSymbols = getChristoffelSymbols(position);
	vec3 velocityDerivative = vec3(0.0);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				velocityDerivative[i] -= christoffelSymbols[i][j][k] * velocity[j] * velocity[k];
			}
		}
	}

	return mat2x3(positionDerivative, velocityDerivative);
}

vec3 sampleUpperBackground(vec3 direction) {
	return Texel(upperSkybox, direction).rgb;
	// return direction * 0.5 + 0.5;
}

vec3 sampleLowerBackground(vec3 direction) {
	return Texel(lowerSkybox, direction).rgb;
	// return (direction * 0.5 + 0.5) * 0.5;
}

void pixelmain() {
	vec3 initialRayDirectionCartesian = normalize(directionPreNormalise);
	bool curvedMode = startInCurvedMode;

	vec3 rayPosition;
	vec3 rayVelocity;
	bool altCoords;
	bool upperWorld;
	if (curvedMode) {
		rayPosition = initialCameraPosition;
		rayVelocity = raySpeed * (cartesianToSphericalBasis(rayPosition) * initialRayDirectionCartesian);
		altCoords = initialAltCoords;
		// upperWorld doesn't need to be set
	} else {
		rayPosition = initialCameraPosition;
		rayVelocity = raySpeed * initialRayDirectionCartesian;
		// altCoords doesn't need to be set
		upperWorld = startInUpperWorld;
	}
	for (int stepNumber = 0; stepNumber < rayStepCount; stepNumber++) {
		float integrationT = 0.0;
		float integrationTStep = 1.0;

		if (!curvedMode) {
			rayPosition += rayVelocity * integrationTStep;

			// Check for flat->curved switch
			float aDistance = distance(rayPosition, mouthAPosition);
			float bDistance = distance(rayPosition, mouthBPosition);
			vec3 mouthPosition = aDistance < bDistance ? mouthAPosition : mouthBPosition;

			vec3 delta = rayPosition - mouthPosition;
			float r = length(delta);
			float cutoffR = getWormholeCutoffR(false);

			if (abs(r) < cutoffR) {
				// TODO (dont forget to set altcoords and use upperworld)
				discard;
				// curvedMode = true;
				// altCoords = false;
			}
		} else {
			// Check if a coordinate transition is needed between main and alternate spherical coords
			if (
				rayPosition.y < altCoordsProportion * tau / 2.0 ||
				rayPosition.y > (1.0 - altCoordsProportion) * tau / 2.0
			) {
				vec3 fakeRayPosition = vec3(1.0, rayPosition.yz); // For some reason setting the 1.0 to a really high value makes things smoother when switching between coordinate systems every frame (for testing purposes), but doesn't deal with the seam seen when switching coordinate systems normally?

				// float r = rayPosition.x;
				// float rho = sqrt(r * r + wormholeThroatRadius * wormholeThroatRadius)
				// vec3 fakeRayPosition = vec3(rho, rayPosition.yz);

				// vec3 fakeRayPosition = rayPosition; // Breaks on the intersection between the wormhole sphere and surface where theta is at a value where we have to switch coordinates. Because r = 0 there

				// if (altCoords) {
				// 	rayVelocity.yz = spherical2ToSpherical1Velocity(fakeRayPosition, rayVelocity).yz;
				// 	rayPosition.yz = spherical2ToSpherical1Position(fakeRayPosition).yz;
				// 	altCoords = false;
				// } else {
				// 	rayVelocity.yz = spherical1ToSpherical2Velocity(fakeRayPosition, rayVelocity).yz;
				// 	rayPosition.yz = spherical1ToSpherical2Position(fakeRayPosition).yz;
				// 	altCoords = true;
				// }

				rayVelocity.yz = spherical1ToSpherical2Velocity(fakeRayPosition, rayVelocity).yz;
				rayPosition.yz = spherical1ToSpherical2Position(fakeRayPosition).yz;
				altCoords = !altCoords;
			}

			// Move (velocity is parallel transported)

			mat2x3 state = mat2x3(
				rayPosition,
				rayVelocity
			);

			// Euler integration
			// state += integrationStateDeriv(integrationT, state) * integrationTStep;

			// Runge-Kutta 4 integration
			mat2x3 k1 = integrationStateDeriv(integrationT, state) * integrationTStep;
			mat2x3 k2 = integrationStateDeriv(integrationT + integrationTStep / 2.0, state + k1 / 2.0) * integrationTStep;
			mat2x3 k3 = integrationStateDeriv(integrationT + integrationTStep / 2.0, state + k2 / 2.0) * integrationTStep;
			mat2x3 k4 = integrationStateDeriv(integrationT + integrationTStep, state + k3) * integrationTStep;
			state += (k1 + k2 * 2.0 + k3 * 2.0 + k4) / 6.0;

			rayPosition = state[0];
			rayVelocity = state[1];

			// Check for a curved->flat mode switch
			float r = rayPosition.x;
			float cutoffR = getWormholeCutoffR(true);
			if (abs(r) > cutoffR) {
				vec3 newPosition;
				vec3 newVelocity;
				if (altCoords) {
					newPosition = spherical2ToAbsoluteCartesianPosition(rayPosition);
					newVelocity = spherical2ToAbsoluteCartesianTangent(rayPosition, rayVelocity);
				} else {
					newPosition = spherical1ToAbsoluteCartesianPosition(rayPosition);
					newVelocity = spherical1ToAbsoluteCartesianTangent(rayPosition, rayVelocity);
				}

				rayPosition = newPosition;
				rayVelocity = newVelocity;
				curvedMode = false;
				upperWorld = r >= 0.0;
			}
		}
	}

	vec3 finalRayDirectionCartesian;
	if (curvedMode) {
		finalRayDirectionCartesian = normalize(altCoords ?
			spherical2ToAbsoluteCartesianTangent(rayPosition, rayVelocity) :
			spherical1ToAbsoluteCartesianTangent(rayPosition, rayVelocity)
		);
		upperWorld = rayPosition.x >= 0;
	} else {
		finalRayDirectionCartesian = normalize(rayVelocity);
		// Retain prior upperWorld variable
	}
	vec3 backgroundColour;
	if (upperWorld) {
		backgroundColour = sampleUpperBackground(finalRayDirectionCartesian);
	} else {
		backgroundColour = sampleLowerBackground(finalRayDirectionCartesian);
	}
	outColour = vec4(backgroundColour, 1.0);
}

#endif

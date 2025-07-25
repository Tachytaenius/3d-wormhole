const float tau = 6.28318530717958647692528676655900576839433879875021164194988918461563281257241799725606965068423413; // Thanks OEIS

uniform vec3 initialCameraForward;
uniform vec3 initialCameraUp;
uniform vec3 initialCameraRight;

#ifdef VERTEX

uniform float cameraHorizontalDirectionExtent;
uniform float cameraVerticalDirectionExtent;

in layout(location = 0) vec4 VertexPosition;

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

uniform vec3 initialCameraPosition;

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

vec3 sampleUpperBackground(vec3 direction) {
	return direction * 0.5 + 0.5;
}

vec3 sampleLowerBackground(vec3 direction) {
	return (direction * 0.5 + 0.5) * 0.5;
}

void pixelmain() {
	vec3 initialRayDirectionCartesian = normalize(directionPreNormalise);
	vec3 rayPosition = initialCameraPosition;
	vec3 rayVelocity = raySpeed * (cartesianToSphericalBasis(rayPosition) * initialRayDirectionCartesian);
	for (int stepNumber = 0; stepNumber < rayStepCount; stepNumber++) {
		// Move (velocity is parallel transported)
		// Euler integration
		mat3[3] christoffelSymbols = getChristoffelSymbols(rayPosition);
		vec3 newPosition = rayPosition + rayVelocity;
		vec3 newVelocity = rayVelocity;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				for (int k = 0; k < 3; k++) {
					newVelocity[i] -= christoffelSymbols[i][j][k] * rayVelocity[j] * rayVelocity[k];
				}
			}
		}
		rayPosition = newPosition;
		rayVelocity = newVelocity;
	}

	vec3 finalRayDirectionCartesian = normalize(sphericalToCartesianBasis(rayPosition) * rayVelocity);
	vec3 backgroundColour;
	if (rayPosition.x > 0) {
		backgroundColour = sampleUpperBackground(finalRayDirectionCartesian);
	} else {
		backgroundColour = sampleLowerBackground(finalRayDirectionCartesian);
	}
	outColour = vec4(backgroundColour, 1.0);
}

#endif

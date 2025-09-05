package com.csc205.project1;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.params.provider.CsvSource;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Comprehensive unit tests for the Point3D class.
 * 
 * This test suite demonstrates several testing principles:
 * - BOUNDARY VALUE TESTING: Testing edge cases and limits
 * - EQUIVALENCE CLASS PARTITIONING: Testing representative values from different input categories
 * - ERROR PATH TESTING: Validating exception handling
 * - STATE-BASED TESTING: Verifying object state after operations
 * - BEHAVIOR-BASED TESTING: Testing method interactions and side effects
 */
@DisplayName("Point3D Tests")
class Point3DTest {

    private static final double EPSILON = 1e-9; // Tolerance for floating-point comparisons
    private Point3D origin;
    private Point3D unitX;
    private Point3D unitY;
    private Point3D unitZ;
    private Point3D testPoint;

    @BeforeEach
    void setUp() {
        origin = new Point3D(0, 0, 0);
        unitX = new Point3D(1, 0, 0);
        unitY = new Point3D(0, 1, 0);
        unitZ = new Point3D(0, 0, 1);
        testPoint = new Point3D(3, 4, 5);
    }

    @Nested
    @DisplayName("Constructor Tests")
    class ConstructorTests {

        @Test
        @DisplayName("Should create point with valid coordinates")
        void shouldCreatePointWithValidCoordinates() {
            Point3D point = new Point3D(1.5, -2.7, 3.14);
            
            assertEquals(1.5, point.getX(), EPSILON);
            assertEquals(-2.7, point.getY(), EPSILON);
            assertEquals(3.14, point.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should create point with zero coordinates")
        void shouldCreatePointWithZeroCoordinates() {
            Point3D point = new Point3D(0, 0, 0);
            
            assertEquals(0, point.getX(), EPSILON);
            assertEquals(0, point.getY(), EPSILON);
            assertEquals(0, point.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should create point with very large coordinates")
        void shouldCreatePointWithVeryLargeCoordinates() {
            double large = 1e100;
            Point3D point = new Point3D(large, -large, large);
            
            assertEquals(large, point.getX(), EPSILON);
            assertEquals(-large, point.getY(), EPSILON);
            assertEquals(large, point.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should create point with very small coordinates")
        void shouldCreatePointWithVerySmallCoordinates() {
            double small = 1e-100;
            Point3D point = new Point3D(small, -small, small);
            
            assertEquals(small, point.getX(), EPSILON);
            assertEquals(-small, point.getY(), EPSILON);
            assertEquals(small, point.getZ(), EPSILON);
        }

        @ParameterizedTest
        @DisplayName("Should reject NaN coordinates")
        @CsvSource({
            "NaN, 0, 0",
            "0, NaN, 0", 
            "0, 0, NaN",
            "NaN, NaN, NaN"
        })
        void shouldRejectNaNCoordinates(String xStr, String yStr, String zStr) {
            double x = "NaN".equals(xStr) ? Double.NaN : Double.parseDouble(xStr);
            double y = "NaN".equals(yStr) ? Double.NaN : Double.parseDouble(yStr);
            double z = "NaN".equals(zStr) ? Double.NaN : Double.parseDouble(zStr);
            
            IllegalArgumentException exception = assertThrows(
                IllegalArgumentException.class,
                () -> new Point3D(x, y, z)
            );
            
            assertEquals("Point coordinates cannot be NaN", exception.getMessage());
        }

        @ParameterizedTest
        @DisplayName("Should reject infinite coordinates")
        @ValueSource(doubles = {Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY})
        void shouldRejectInfiniteCoordinates(double infiniteValue) {
            IllegalArgumentException exception1 = assertThrows(
                IllegalArgumentException.class,
                () -> new Point3D(infiniteValue, 0, 0)
            );
            
            IllegalArgumentException exception2 = assertThrows(
                IllegalArgumentException.class,
                () -> new Point3D(0, infiniteValue, 0)
            );
            
            IllegalArgumentException exception3 = assertThrows(
                IllegalArgumentException.class,
                () -> new Point3D(0, 0, infiniteValue)
            );
            
            assertEquals("Point coordinates cannot be infinite", exception1.getMessage());
            assertEquals("Point coordinates cannot be infinite", exception2.getMessage());
            assertEquals("Point coordinates cannot be infinite", exception3.getMessage());
        }
    }

    @Nested
    @DisplayName("Static Factory Method Tests")
    class StaticFactoryTests {

        @Test
        @DisplayName("Should verify predefined constants")
        void shouldVerifyPredefinedConstants() {
            assertEquals(0, Point3D.ORIGIN.getX(), EPSILON);
            assertEquals(0, Point3D.ORIGIN.getY(), EPSILON);
            assertEquals(0, Point3D.ORIGIN.getZ(), EPSILON);
            
            assertEquals(1, Point3D.UNIT_X.getX(), EPSILON);
            assertEquals(0, Point3D.UNIT_X.getY(), EPSILON);
            assertEquals(0, Point3D.UNIT_X.getZ(), EPSILON);
            
            assertEquals(0, Point3D.UNIT_Y.getX(), EPSILON);
            assertEquals(1, Point3D.UNIT_Y.getY(), EPSILON);
            assertEquals(0, Point3D.UNIT_Y.getZ(), EPSILON);
            
            assertEquals(0, Point3D.UNIT_Z.getX(), EPSILON);
            assertEquals(0, Point3D.UNIT_Z.getY(), EPSILON);
            assertEquals(1, Point3D.UNIT_Z.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should create point from spherical coordinates")
        void shouldCreatePointFromSphericalCoordinates() {
            double radius = 5.0;
            double theta = Math.PI / 4; // 45 degrees
            double phi = Math.PI / 3;   // 60 degrees
            
            Point3D point = Point3D.fromSpherical(radius, theta, phi);
            
            double expectedX = radius * Math.sin(phi) * Math.cos(theta);
            double expectedY = radius * Math.sin(phi) * Math.sin(theta);
            double expectedZ = radius * Math.cos(phi);
            
            assertEquals(expectedX, point.getX(), EPSILON);
            assertEquals(expectedY, point.getY(), EPSILON);
            assertEquals(expectedZ, point.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should create origin from zero radius spherical coordinates")
        void shouldCreateOriginFromZeroRadiusSphericalCoordinates() {
            Point3D point = Point3D.fromSpherical(0, Math.PI, Math.PI / 2);
            
            assertEquals(0, point.getX(), EPSILON);
            assertEquals(0, point.getY(), EPSILON);
            assertEquals(0, point.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should reject negative radius in spherical coordinates")
        void shouldRejectNegativeRadiusInSphericalCoordinates() {
            IllegalArgumentException exception = assertThrows(
                IllegalArgumentException.class,
                () -> Point3D.fromSpherical(-1.0, 0, 0)
            );
            
            assertEquals("Radius must be non-negative", exception.getMessage());
        }

        @Test
        @DisplayName("Should handle extreme spherical coordinate values")
        void shouldHandleExtremeSphericalCoordinateValues() {
            // Test with very large radius
            Point3D largeRadius = Point3D.fromSpherical(1e10, 0, Math.PI / 2);
            assertTrue(Math.abs(largeRadius.getX()) > 1e9);
            
            // Test with very small radius  
            Point3D smallRadius = Point3D.fromSpherical(1e-10, Math.PI, Math.PI / 2);
            assertTrue(Math.abs(smallRadius.getX()) < 1e-9);
        }
    }

    @Nested
    @DisplayName("Distance Calculation Tests")
    class DistanceCalculationTests {

        @Test
        @DisplayName("Should calculate distance between two points")
        void shouldCalculateDistanceBetweenTwoPoints() {
            Point3D point1 = new Point3D(0, 0, 0);
            Point3D point2 = new Point3D(3, 4, 0);
            
            double distance = point1.distanceTo(point2);
            assertEquals(5.0, distance, EPSILON);
        }

        @Test
        @DisplayName("Should calculate distance in 3D space")
        void shouldCalculateDistanceIn3DSpace() {
            Point3D point1 = new Point3D(1, 2, 3);
            Point3D point2 = new Point3D(4, 6, 8);
            
            double expectedDistance = Math.sqrt(9 + 16 + 25); // sqrt(50)
            double actualDistance = point1.distanceTo(point2);
            
            assertEquals(expectedDistance, actualDistance, EPSILON);
        }

        @Test
        @DisplayName("Should return zero distance for identical points")
        void shouldReturnZeroDistanceForIdenticalPoints() {
            double distance = testPoint.distanceTo(testPoint);
            assertEquals(0.0, distance, EPSILON);
        }

        @Test
        @DisplayName("Should calculate symmetric distance")
        void shouldCalculateSymmetricDistance() {
            double distance1 = testPoint.distanceTo(origin);
            double distance2 = origin.distanceTo(testPoint);
            
            assertEquals(distance1, distance2, EPSILON);
        }

        @Test
        @DisplayName("Should throw exception for null point in distance calculation")
        void shouldThrowExceptionForNullPointInDistanceCalculation() {
            IllegalArgumentException exception = assertThrows(
                IllegalArgumentException.class,
                () -> testPoint.distanceTo(null)
            );
            
            assertEquals("Other point cannot be null", exception.getMessage());
        }

        @Test
        @DisplayName("Should calculate squared distance correctly")
        void shouldCalculateSquaredDistanceCorrectly() {
            Point3D point1 = new Point3D(1, 2, 3);
            Point3D point2 = new Point3D(4, 6, 8);
            
            double squaredDistance = point1.distanceSquaredTo(point2);
            double distance = point1.distanceTo(point2);
            
            assertEquals(distance * distance, squaredDistance, EPSILON);
        }

        @Test
        @DisplayName("Should handle very large distances")
        void shouldHandleVeryLargeDistances() {
            Point3D point1 = new Point3D(0, 0, 0);
            Point3D point2 = new Point3D(1e50, 1e50, 1e50);
            
            double distance = point1.distanceTo(point2);
            assertTrue(distance > 1e50);
            assertFalse(Double.isInfinite(distance));
        }
    }

    @Nested
    @DisplayName("Rotation Tests")
    class RotationTests {

        @Test
        @DisplayName("Should rotate around X-axis correctly")
        void shouldRotateAroundXAxisCorrectly() {
            Point3D point = new Point3D(1, 1, 0);
            Point3D rotated = point.rotateX(Math.PI / 2); // 90 degrees
            
            assertEquals(1, rotated.getX(), EPSILON);
            assertEquals(0, rotated.getY(), EPSILON);
            assertEquals(1, rotated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should rotate around Y-axis correctly")
        void shouldRotateAroundYAxisCorrectly() {
            Point3D point = new Point3D(1, 1, 0);
            Point3D rotated = point.rotateY(Math.PI / 2); // 90 degrees
            
            assertEquals(0, rotated.getX(), EPSILON);
            assertEquals(1, rotated.getY(), EPSILON);
            assertEquals(-1, rotated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should rotate around Z-axis correctly")
        void shouldRotateAroundZAxisCorrectly() {
            Point3D point = new Point3D(1, 0, 1);
            Point3D rotated = point.rotateZ(Math.PI / 2); // 90 degrees
            
            assertEquals(0, rotated.getX(), EPSILON);
            assertEquals(1, rotated.getY(), EPSILON);
            assertEquals(1, rotated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should preserve distance from origin during rotation")
        void shouldPreserveDistanceFromOriginDuringRotation() {
            Point3D original = new Point3D(3, 4, 5);
            double originalDistance = original.distanceTo(origin);
            
            Point3D rotatedX = original.rotateX(Math.PI / 3);
            Point3D rotatedY = original.rotateY(Math.PI / 4);
            Point3D rotatedZ = original.rotateZ(Math.PI / 6);
            
            assertEquals(originalDistance, rotatedX.distanceTo(origin), EPSILON);
            assertEquals(originalDistance, rotatedY.distanceTo(origin), EPSILON);
            assertEquals(originalDistance, rotatedZ.distanceTo(origin), EPSILON);
        }

        @Test
        @DisplayName("Should return same point for zero rotation")
        void shouldReturnSamePointForZeroRotation() {
            Point3D rotatedX = testPoint.rotateX(0);
            Point3D rotatedY = testPoint.rotateY(0);
            Point3D rotatedZ = testPoint.rotateZ(0);
            
            assertEquals(testPoint.getX(), rotatedX.getX(), EPSILON);
            assertEquals(testPoint.getY(), rotatedX.getY(), EPSILON);
            assertEquals(testPoint.getZ(), rotatedX.getZ(), EPSILON);
            
            assertEquals(testPoint, rotatedY);
            assertEquals(testPoint, rotatedZ);
        }

        @Test
        @DisplayName("Should handle full rotation (2π)")
        void shouldHandleFullRotation() {
            Point3D rotated = testPoint.rotateX(2 * Math.PI);
            
            assertEquals(testPoint.getX(), rotated.getX(), EPSILON);
            assertEquals(testPoint.getY(), rotated.getY(), EPSILON);
            assertEquals(testPoint.getZ(), rotated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should return original point for invalid rotation angle")
        void shouldReturnOriginalPointForInvalidRotationAngle() {
            Point3D rotatedNaN = testPoint.rotateX(Double.NaN);
            Point3D rotatedInf = testPoint.rotateY(Double.POSITIVE_INFINITY);
            
            assertEquals(testPoint, rotatedNaN);
            assertEquals(testPoint, rotatedInf);
        }
    }

    @Nested
    @DisplayName("Translation Tests")
    class TranslationTests {

        @Test
        @DisplayName("Should translate point correctly")
        void shouldTranslatePointCorrectly() {
            Point3D translated = testPoint.translate(1, -2, 3);
            
            assertEquals(4, translated.getX(), EPSILON);
            assertEquals(2, translated.getY(), EPSILON);
            assertEquals(8, translated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should return same point for zero translation")
        void shouldReturnSamePointForZeroTranslation() {
            Point3D translated = testPoint.translate(0, 0, 0);
            assertEquals(testPoint, translated);
        }

        @Test
        @DisplayName("Should handle negative translation values")
        void shouldHandleNegativeTranslationValues() {
            Point3D translated = origin.translate(-1, -2, -3);
            
            assertEquals(-1, translated.getX(), EPSILON);
            assertEquals(-2, translated.getY(), EPSILON);
            assertEquals(-3, translated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should handle very large translation values")
        void shouldHandleVeryLargeTranslationValues() {
            double large = 1e100;
            Point3D translated = origin.translate(large, -large, large);
            
            assertEquals(large, translated.getX(), EPSILON);
            assertEquals(-large, translated.getY(), EPSILON);
            assertEquals(large, translated.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should return original point for NaN translation")
        void shouldReturnOriginalPointForNaNTranslation() {
            Point3D translatedNaN = testPoint.translate(Double.NaN, 1, 2);
            Point3D translatedNaN2 = testPoint.translate(1, Double.NaN, 2);
            Point3D translatedNaN3 = testPoint.translate(1, 2, Double.NaN);
            
            assertEquals(testPoint, translatedNaN);
            assertEquals(testPoint, translatedNaN2);
            assertEquals(testPoint, translatedNaN3);
        }
    }

    @Nested
    @DisplayName("Vector Operations Tests")
    class VectorOperationTests {

        @Test
        @DisplayName("Should calculate dot product correctly")
        void shouldCalculateDotProductCorrectly() {
            Point3D point1 = new Point3D(1, 2, 3);
            Point3D point2 = new Point3D(4, 5, 6);
            
            double dotProduct = point1.dotProduct(point2);
            double expected = 1*4 + 2*5 + 3*6; // 4 + 10 + 18 = 32
            
            assertEquals(expected, dotProduct, EPSILON);
        }

        @Test
        @DisplayName("Should return zero for perpendicular vectors dot product")
        void shouldReturnZeroForPerpendicularVectorsDotProduct() {
            Point3D perpendicular1 = new Point3D(1, 0, 0);
            Point3D perpendicular2 = new Point3D(0, 1, 0);
            
            double dotProduct = perpendicular1.dotProduct(perpendicular2);
            assertEquals(0, dotProduct, EPSILON);
        }

        @Test
        @DisplayName("Should calculate cross product correctly")
        void shouldCalculateCrossProductCorrectly() {
            Point3D point1 = new Point3D(1, 0, 0);
            Point3D point2 = new Point3D(0, 1, 0);
            
            Point3D crossProduct = point1.crossProduct(point2);
            
            assertEquals(0, crossProduct.getX(), EPSILON);
            assertEquals(0, crossProduct.getY(), EPSILON);
            assertEquals(1, crossProduct.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should return zero vector for parallel vectors cross product")
        void shouldReturnZeroVectorForParallelVectorsCrossProduct() {
            Point3D vector1 = new Point3D(2, 4, 6);
            Point3D vector2 = new Point3D(1, 2, 3); // parallel to vector1
            
            Point3D crossProduct = vector1.crossProduct(vector2);
            
            assertEquals(0, crossProduct.getX(), EPSILON);
            assertEquals(0, crossProduct.getY(), EPSILON);
            assertEquals(0, crossProduct.getZ(), EPSILON);
        }

        @Test
        @DisplayName("Should throw exception for null point in vector operations")
        void shouldThrowExceptionForNullPointInVectorOperations() {
            assertThrows(IllegalArgumentException.class, () -> testPoint.dotProduct(null));
            assertThrows(IllegalArgumentException.class, () -> testPoint.crossProduct(null));
        }

        @Test
        @DisplayName("Should calculate magnitude correctly")
        void shouldCalculateMagnitudeCorrectly() {
            Point3D point = new Point3D(3, 4, 0);
            double magnitude = point.magnitude();
            assertEquals(5.0, magnitude, EPSILON);
        }

        @Test
        @DisplayName("Should return zero magnitude for origin")
        void shouldReturnZeroMagnitudeForOrigin() {
            double magnitude = origin.magnitude();
            assertEquals(0.0, magnitude, EPSILON);
        }

        @Test
        @DisplayName("Should normalize vector correctly")
        void shouldNormalizeVectorCorrectly() {
            Point3D point = new Point3D(3, 4, 0);
            Point3D normalized = point.normalize();
            
            assertEquals(0.6, normalized.getX(), EPSILON);
            assertEquals(0.8, normalized.getY(), EPSILON);
            assertEquals(0.0, normalized.getZ(), EPSILON);
            
            // Normalized vector should have magnitude 1
            assertEquals(1.0, normalized.magnitude(), EPSILON);
        }

        @Test
        @DisplayName("Should throw exception when normalizing zero vector")
        void shouldThrowExceptionWhenNormalizingZeroVector() {
            IllegalStateException exception = assertThrows(
                IllegalStateException.class,
                () -> origin.normalize()
            );
            
            assertEquals("Cannot normalize zero vector", exception.getMessage());
        }
    }

    @Nested
    @DisplayName("Equality and Hash Code Tests")
    class EqualityAndHashCodeTests {

        @Test
        @DisplayName("Should be equal to itself")
        void shouldBeEqualToItself() {
            assertEquals(testPoint, testPoint);
        }

        @Test
        @DisplayName("Should be equal to point with same coordinates")
        void shouldBeEqualToPointWithSameCoordinates() {
            Point3D samePoint = new Point3D(3, 4, 5);
            assertEquals(testPoint, samePoint);
        }

        @Test
        @DisplayName("Should not be equal to point with different coordinates")
        void shouldNotBeEqualToPointWithDifferentCoordinates() {
            Point3D differentPoint = new Point3D(3, 4, 6);
            assertNotEquals(testPoint, differentPoint);
        }

        @Test
        @DisplayName("Should not be equal to null")
        void shouldNotBeEqualToNull() {
            assertNotEquals(testPoint, null);
        }

        @Test
        @DisplayName("Should not be equal to different class")
        void shouldNotBeEqualToDifferentClass() {
            assertNotEquals(testPoint, "Not a Point3D");
        }

        @Test
        @DisplayName("Should have same hash code for equal points")
        void shouldHaveSameHashCodeForEqualPoints() {
            Point3D samePoint = new Point3D(3, 4, 5);
            assertEquals(testPoint.hashCode(), samePoint.hashCode());
        }

        @Test
        @DisplayName("Should handle equality with very small differences")
        void shouldHandleEqualityWithVerySmallDifferences() {
            Point3D point1 = new Point3D(1.0, 2.0, 3.0);
            Point3D point2 = new Point3D(1.0000000000001, 2.0, 3.0);
            
            // These should not be equal due to floating-point precision
            assertNotEquals(point1, point2);
        }
    }

    @Nested
    @DisplayName("String Representation Tests")
    class StringRepresentationTests {

        @Test
        @DisplayName("Should format string representation correctly")
        void shouldFormatStringRepresentationCorrectly() {
            String str = testPoint.toString();
            assertEquals("Point3D(3.000, 4.000, 5.000)", str);
        }

        @Test
        @DisplayName("Should format origin correctly")
        void shouldFormatOriginCorrectly() {
            String str = origin.toString();
            assertEquals("Point3D(0.000, 0.000, 0.000)", str);
        }

        @Test
        @DisplayName("Should format negative coordinates correctly")
        void shouldFormatNegativeCoordinatesCorrectly() {
            Point3D negativePoint = new Point3D(-1.5, -2.7, -3.14);
            String str = negativePoint.toString();
            assertEquals("Point3D(-1.500, -2.700, -3.140)", str);
        }
    }

    @Nested
    @DisplayName("Method Chaining Tests")
    class MethodChainingTests {

        @Test
        @DisplayName("Should support method chaining for transformations")
        void shouldSupportMethodChainingForTransformations() {
            Point3D result = testPoint
                .translate(1, 1, 1)
                .rotateX(Math.PI / 2)
                .rotateY(Math.PI / 4)
                .translate(-1, -1, -1);
            
            assertNotNull(result);
            // Verify that the result is a valid Point3D
            assertFalse(Double.isNaN(result.getX()));
            assertFalse(Double.isNaN(result.getY()));
            assertFalse(Double.isNaN(result.getZ()));
        }

        @Test
        @DisplayName("Should maintain immutability during method chaining")
        void shouldMaintainImmutabilityDuringMethodChaining() {
            Point3D original = new Point3D(1, 2, 3);
            Point3D transformed = original.translate(1, 1, 1).rotateZ(Math.PI);
            
            // Original should remain unchanged
            assertEquals(1, original.getX(), EPSILON);
            assertEquals(2, original.getY(), EPSILON);
            assertEquals(3, original.getZ(), EPSILON);
            
            // Transformed should be different
            assertNotEquals(original, transformed);
        }
    }

    @Nested
    @DisplayName("Edge Case Tests")
    class EdgeCaseTests {

        @Test
        @DisplayName("Should handle operations with very small numbers")
        void shouldHandleOperationsWithVerySmallNumbers() {
            Point3D tiny = new Point3D(1e-100, 1e-100, 1e-100);
            Point3D result = tiny.normalize();
            
            assertEquals(1.0, result.magnitude(), 1e-10);
        }

        @Test
        @DisplayName("Should handle operations with very large numbers")  
        void shouldHandleOperationsWithVeryLargeNumbers() {
            Point3D huge = new Point3D(1e50, 1e50, 1e50);
            double magnitude = huge.magnitude();
            
            assertTrue(magnitude > 1e50);
            assertFalse(Double.isInfinite(magnitude));
        }

        @Test
        @DisplayName("Should handle mixed large and small coordinates")
        void shouldHandleMixedLargeAndSmallCoordinates() {
            Point3D mixed = new Point3D(1e-100, 1e100, 1e-50);
            double magnitude = mixed.magnitude();
            
            // Should be dominated by the large coordinate
            assertTrue(magnitude > 1e99);
        }

        @Test
        @DisplayName("Should maintain precision in complex operations")
        void shouldMaintainPrecisionInComplexOperations() {
            Point3D point = new Point3D(1, 0, 0);
            
            // Rotate 4 times by π/2 around Z-axis should return to original
            Point3D result = point
                .rotateZ(Math.PI / 2)
                .rotateZ(Math.PI / 2)
                .rotateZ(Math.PI / 2)
                .rotateZ(Math.PI / 2);
            
            assertEquals(1, result.getX(), 1e-10);
            assertEquals(0, result.getY(), 1e-10);
            assertEquals(0, result.getZ(), 1e-10);
        }
    }
}
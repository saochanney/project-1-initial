package com.csc205.project1;

import java.util.Objects;
import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a point in 3D space with x, y, and z coordinates.
 * 
 * This class demonstrates several fundamental object-oriented design patterns:
 * 
 * 1. VALUE OBJECT PATTERN: Point3D is immutable, ensuring thread safety and predictable behavior.
 *    All operations return new instances rather than modifying the existing object.
 * 
 * 2. BUILDER-STYLE FACTORY METHODS: Static factory methods provide clear, readable ways to create points.
 * 
 * 3. FLUENT INTERFACE: Method chaining allows for readable transformation sequences.
 * 
 * 4. TEMPLATE METHOD PATTERN: Common mathematical operations are structured consistently.
 * 
 * 5. DEFENSIVE PROGRAMMING: Input validation and error handling throughout.
 * 
 * Data Structure Principles Demonstrated:
 * - Encapsulation: Coordinates are private with controlled access
 * - Immutability: Thread-safe design prevents concurrent modification issues
 * - Memory efficiency: Primitive doubles used for coordinates
 * 
 * Algorithm Principles Demonstrated:
 * - Vector mathematics: Distance, dot product, cross product calculations
 * - Linear algebra: Rotation matrices and transformations
 * - Geometric algorithms: Point-to-point and point-to-plane operations
 * - Numerical stability: Careful handling of floating-point operations
 * 
 * @author Spring Framework Style Documentation
 * @version 1.0
 */
public final class Point3D {
    
    private static final Logger logger = Logger.getLogger(Point3D.class.getName());
    
    // Immutable coordinates - demonstrates VALUE OBJECT pattern
    private final double x;
    private final double y;
    private final double z;
    
    // Common point constants - demonstrates FLYWEIGHT pattern for frequently used objects
    public static final Point3D ORIGIN = new Point3D(0, 0, 0);
    public static final Point3D UNIT_X = new Point3D(1, 0, 0);
    public static final Point3D UNIT_Y = new Point3D(0, 1, 0);
    public static final Point3D UNIT_Z = new Point3D(0, 0, 1);
    
    /**
     * Constructs a new Point3D with the specified coordinates.
     * 
     * This constructor follows the Spring Framework philosophy of explicit parameter validation
     * and clear error messaging. It performs input validation to ensure numerical stability
     * and logs the creation of new points for debugging purposes.
     * 
     * The constructor demonstrates the DEFENSIVE PROGRAMMING pattern by validating inputs
     * and the IMMUTABLE OBJECT pattern by making all fields final.
     * 
     * @param x the x-coordinate of the point
     * @param y the y-coordinate of the point  
     * @param z the z-coordinate of the point
     * @throws IllegalArgumentException if any coordinate is NaN or infinite
     */
    public Point3D(double x, double y, double z) {
        // Input validation - demonstrates defensive programming
        if (Double.isNaN(x) || Double.isNaN(y) || Double.isNaN(z)) {
            logger.log(Level.SEVERE, "Attempted to create Point3D with NaN coordinates: x={0}, y={1}, z={2}", 
                      new Object[]{x, y, z});
            throw new IllegalArgumentException("Point coordinates cannot be NaN");
        }
        
        if (Double.isInfinite(x) || Double.isInfinite(y) || Double.isInfinite(z)) {
            logger.log(Level.SEVERE, "Attempted to create Point3D with infinite coordinates: x={0}, y={1}, z={2}", 
                      new Object[]{x, y, z});
            throw new IllegalArgumentException("Point coordinates cannot be infinite");
        }
        
        this.x = x;
        this.y = y;
        this.z = z;
        
        logger.log(Level.INFO, "Created new Point3D: ({0}, {1}, {2})", new Object[]{x, y, z});
    }
    
    /**
     * Static factory method to create a Point3D from spherical coordinates.
     * 
     * This method demonstrates the FACTORY METHOD pattern, providing a clear and intuitive
     * way to create points using spherical coordinates (radius, theta, phi). This approach
     * follows Spring's convention of using static factory methods with descriptive names
     * rather than overloaded constructors.
     * 
     * The conversion uses standard mathematical formulas:
     * - x = radius * sin(phi) * cos(theta)
     * - y = radius * sin(phi) * sin(theta)  
     * - z = radius * cos(phi)
     * 
     * @param radius the distance from origin (must be non-negative)
     * @param theta the azimuthal angle in radians
     * @param phi the polar angle in radians
     * @return a new Point3D representing the spherical coordinates
     * @throws IllegalArgumentException if radius is negative
     */
    public static Point3D fromSpherical(double radius, double theta, double phi) {
        if (radius < 0) {
            logger.log(Level.SEVERE, "Invalid spherical coordinates: negative radius {0}", radius);
            throw new IllegalArgumentException("Radius must be non-negative");
        }
        
        double x = radius * Math.sin(phi) * Math.cos(theta);
        double y = radius * Math.sin(phi) * Math.sin(theta);
        double z = radius * Math.cos(phi);
        
        logger.log(Level.INFO, "Created Point3D from spherical coordinates: r={0}, θ={1}, φ={2}", 
                  new Object[]{radius, theta, phi});
        
        return new Point3D(x, y, z);
    }
    
    /**
     * Calculates the Euclidean distance between this point and another point.
     * 
     * This method implements the fundamental distance formula in 3D space:
     * distance = √((x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²)
     * 
     * The implementation demonstrates several algorithmic principles:
     * - Numerical stability through careful ordering of operations
     * - Input validation following the DEFENSIVE PROGRAMMING pattern
     * - Clear error handling and logging
     * 
     * Time complexity: O(1) - constant time operation
     * Space complexity: O(1) - uses only a constant amount of extra space
     * 
     * @param other the other point to calculate distance to (must not be null)
     * @return the Euclidean distance between this point and the other point
     * @throws IllegalArgumentException if the other point is null
     */
    public double distanceTo(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot calculate distance: other point is null");
            throw new IllegalArgumentException("Other point cannot be null");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
        
        logger.log(Level.INFO, "Calculated distance from {0} to {1}: {2}", 
                  new Object[]{this, other, distance});
        
        return distance;
    }
    
    /**
     * Calculates the squared distance between this point and another point.
     * 
     * This method provides an optimized distance calculation when you only need to compare
     * distances (avoiding the expensive square root operation). This optimization is commonly
     * used in algorithms like k-nearest neighbors or spatial partitioning.
     * 
     * This demonstrates the PERFORMANCE OPTIMIZATION pattern where we provide alternative
     * implementations for performance-critical scenarios.
     * 
     * Time complexity: O(1) - constant time, faster than distanceTo()
     * Space complexity: O(1) - constant space
     * 
     * @param other the other point to calculate squared distance to
     * @return the squared Euclidean distance
     * @throws IllegalArgumentException if the other point is null
     */
    public double distanceSquaredTo(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot calculate squared distance: other point is null");
            throw new IllegalArgumentException("Other point cannot be null");
        }
        
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        double dz = this.z - other.z;
        
        double distanceSquared = dx * dx + dy * dy + dz * dz;
        
        logger.log(Level.INFO, "Calculated squared distance from {0} to {1}: {2}", 
                  new Object[]{this, other, distanceSquared});
        
        return distanceSquared;
    }
    
    /**
     * Rotates this point around the X-axis by the specified angle.
     * 
     * This method implements 3D rotation using rotation matrices, a fundamental concept
     * in computer graphics and linear algebra. The rotation matrix for X-axis rotation is:
     * 
     * [1   0      0    ]
     * [0  cos(θ) -sin(θ)]  
     * [0  sin(θ)  cos(θ)]
     * 
     * This demonstrates the IMMUTABLE OBJECT pattern - instead of modifying this point,
     * we return a new Point3D with the rotated coordinates.
     * 
     * @param angleRadians the rotation angle in radians (positive = counterclockwise when looking from positive X toward origin)
     * @return a new Point3D representing this point rotated around the X-axis
     */
    public Point3D rotateX(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.log(Level.WARNING, "Invalid rotation angle: {0}, returning original point", angleRadians);
            return this;
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newY = this.y * cos - this.z * sin;
        double newZ = this.y * sin + this.z * cos;
        
        logger.log(Level.INFO, "Rotated point {0} around X-axis by {1} radians", 
                  new Object[]{this, angleRadians});
        
        return new Point3D(this.x, newY, newZ);
    }
    
    /**
     * Rotates this point around the Y-axis by the specified angle.
     * 
     * Implements Y-axis rotation using the rotation matrix:
     * [cos(θ)   0  sin(θ)]
     * [0        1  0     ]
     * [-sin(θ)  0  cos(θ)]
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing this point rotated around the Y-axis
     */
    public Point3D rotateY(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.log(Level.WARNING, "Invalid rotation angle: {0}, returning original point", angleRadians);
            return this;
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = this.x * cos + this.z * sin;
        double newZ = -this.x * sin + this.z * cos;
        
        logger.log(Level.INFO, "Rotated point {0} around Y-axis by {1} radians", 
                  new Object[]{this, angleRadians});
        
        return new Point3D(newX, this.y, newZ);
    }
    
    /**
     * Rotates this point around the Z-axis by the specified angle.
     * 
     * Implements Z-axis rotation using the rotation matrix:
     * [cos(θ) -sin(θ)  0]
     * [sin(θ)  cos(θ)  0]
     * [0       0       1]
     * 
     * @param angleRadians the rotation angle in radians
     * @return a new Point3D representing this point rotated around the Z-axis
     */
    public Point3D rotateZ(double angleRadians) {
        if (Double.isNaN(angleRadians) || Double.isInfinite(angleRadians)) {
            logger.log(Level.WARNING, "Invalid rotation angle: {0}, returning original point", angleRadians);
            return this;
        }
        
        double cos = Math.cos(angleRadians);
        double sin = Math.sin(angleRadians);
        
        double newX = this.x * cos - this.y * sin;
        double newY = this.x * sin + this.y * cos;
        
        logger.log(Level.INFO, "Rotated point {0} around Z-axis by {1} radians", 
                  new Object[]{this, angleRadians});
        
        return new Point3D(newX, newY, this.z);
    }
    
    /**
     * Translates this point by the specified vector amounts.
     * 
     * This method implements vector addition, a fundamental operation in linear algebra.
     * Translation is one of the basic geometric transformations along with rotation and scaling.
     * 
     * The method demonstrates the IMMUTABLE OBJECT pattern and FLUENT INTERFACE design,
     * allowing for method chaining: point.translate(1, 0, 0).rotateX(Math.PI/4)
     * 
     * @param dx the amount to translate along the X-axis
     * @param dy the amount to translate along the Y-axis  
     * @param dz the amount to translate along the Z-axis
     * @return a new Point3D representing this point translated by the specified amounts
     */
    public Point3D translate(double dx, double dy, double dz) {
        if (Double.isNaN(dx) || Double.isNaN(dy) || Double.isNaN(dz)) {
            logger.log(Level.WARNING, "Invalid translation values contain NaN: dx={0}, dy={1}, dz={2}", 
                      new Object[]{dx, dy, dz});
            return this;
        }
        
        logger.log(Level.INFO, "Translating point {0} by vector ({1}, {2}, {3})", 
                  new Object[]{this, dx, dy, dz});
        
        return new Point3D(this.x + dx, this.y + dy, this.z + dz);
    }
    
    /**
     * Calculates the dot product of this point (treated as a vector) with another point.
     * 
     * The dot product is a fundamental operation in linear algebra with many applications:
     * - Calculating angles between vectors: cos(θ) = (a·b) / (|a||b|)
     * - Determining orthogonality: vectors are perpendicular if dot product = 0
     * - Projecting one vector onto another
     * 
     * Formula: a·b = ax*bx + ay*by + az*bz
     * 
     * This method demonstrates treating points as position vectors from the origin,
     * a common technique in computer graphics and physics simulations.
     * 
     * @param other the other point to calculate dot product with
     * @return the dot product of the two vectors
     * @throws IllegalArgumentException if the other point is null
     */
    public double dotProduct(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot calculate dot product: other point is null");
            throw new IllegalArgumentException("Other point cannot be null");
        }
        
        double result = this.x * other.x + this.y * other.y + this.z * other.z;
        
        logger.log(Level.INFO, "Calculated dot product of {0} and {1}: {2}", 
                  new Object[]{this, other, result});
        
        return result;
    }
    
    /**
     * Calculates the cross product of this point (as vector) with another point.
     * 
     * The cross product produces a vector perpendicular to both input vectors.
     * It's essential for calculating surface normals, determining rotation directions,
     * and many physics calculations.
     * 
     * Formula: a × b = (ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx)
     * 
     * Properties of cross product:
     * - Result is perpendicular to both input vectors
     * - Magnitude equals area of parallelogram formed by the vectors
     * - Direction follows right-hand rule
     * 
     * @param other the other point to calculate cross product with
     * @return a new Point3D representing the cross product vector
     * @throws IllegalArgumentException if the other point is null
     */
    public Point3D crossProduct(Point3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot calculate cross product: other point is null");
            throw new IllegalArgumentException("Other point cannot be null");
        }
        
        double newX = this.y * other.z - this.z * other.y;
        double newY = this.z * other.x - this.x * other.z;
        double newZ = this.x * other.y - this.y * other.x;
        
        Point3D result = new Point3D(newX, newY, newZ);
        
        logger.log(Level.INFO, "Calculated cross product of {0} and {1}: {2}", 
                  new Object[]{this, other, result});
        
        return result;
    }
    
    /**
     * Calculates the magnitude (length) of this point when treated as a vector from origin.
     * 
     * This method computes the Euclidean norm of the position vector:
     * magnitude = √(x² + y² + z²)
     * 
     * The magnitude is used in vector normalization, distance calculations,
     * and determining vector lengths for physics simulations.
     * 
     * @return the magnitude of this point as a position vector
     */
    public double magnitude() {
        double mag = Math.sqrt(x * x + y * y + z * z);
        
        logger.log(Level.INFO, "Calculated magnitude of {0}: {1}", new Object[]{this, mag});
        
        return mag;
    }
    
    /**
     * Returns a normalized version of this point (unit vector in same direction).
     * 
     * Normalization converts a vector to unit length while preserving its direction.
     * This is crucial for many algorithms including lighting calculations,
     * collision detection, and camera controls.
     * 
     * Formula: normalized = vector / |vector|
     * 
     * @return a new Point3D representing the normalized vector
     * @throws IllegalStateException if this point is at the origin (cannot normalize zero vector)
     */
    public Point3D normalize() {
        double mag = magnitude();
        
        if (mag == 0.0) {
            logger.log(Level.SEVERE, "Cannot normalize zero vector: {0}", this);
            throw new IllegalStateException("Cannot normalize zero vector");
        }
        
        Point3D normalized = new Point3D(x / mag, y / mag, z / mag);
        
        logger.log(Level.INFO, "Normalized vector {0} to {1}", new Object[]{this, normalized});
        
        return normalized;
    }
    
    // Accessor methods - demonstrate ENCAPSULATION principle
    
    /**
     * Returns the X coordinate of this point.
     * @return the X coordinate
     */
    public double getX() { return x; }
    
    /**
     * Returns the Y coordinate of this point.
     * @return the Y coordinate  
     */
    public double getY() { return y; }
    
    /**
     * Returns the Z coordinate of this point.
     * @return the Z coordinate
     */
    public double getZ() { return z; }
    
    /**
     * Indicates whether some other object is "equal to" this one.
     * 
     * This implementation follows the contract specified by Object.equals():
     * - Reflexive: x.equals(x) returns true
     * - Symmetric: x.equals(y) returns true if and only if y.equals(x) returns true
     * - Transitive: if x.equals(y) and y.equals(z), then x.equals(z)
     * - Consistent: multiple invocations return the same result
     * - null handling: x.equals(null) returns false
     * 
     * Uses Double.compare() to handle NaN and infinity cases correctly.
     * 
     * @param obj the reference object with which to compare
     * @return true if this object is equal to the obj argument; false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Point3D point3D = (Point3D) obj;
        return Double.compare(point3D.x, x) == 0 &&
               Double.compare(point3D.y, y) == 0 &&
               Double.compare(point3D.z, z) == 0;
    }
    
    /**
     * Returns a hash code value for the object.
     * 
     * This implementation satisfies the general contract of hashCode():
     * - If two objects are equal according to equals(), they must have the same hash code
     * - Hash code should remain consistent during object lifetime
     * - Uses Objects.hash() for robust handling of double values
     * 
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return Objects.hash(x, y, z);
    }
    
    /**
     * Returns a string representation of this point.
     * 
     * The format is "Point3D(x, y, z)" which is clear and follows common conventions
     * for mathematical objects. This aids in debugging and logging.
     * 
     * @return a string representation of the object
     */
    @Override
    public String toString() {
        return String.format("Point3D(%.3f, %.3f, %.3f)", x, y, z);
    }
}
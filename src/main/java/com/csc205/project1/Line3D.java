package com.csc205.project1;

import java.util.Objects;
import java.util.Optional;
import java.util.logging.Logger;
import java.util.logging.Level;

/**
 * Represents a line segment in 3D space defined by two endpoints.
 * 
 * This class demonstrates several advanced object-oriented design patterns and geometric algorithms:
 * 
 * 1. COMPOSITION PATTERN: Line3D is composed of two Point3D objects, demonstrating "has-a" relationships
 *    rather than inheritance. This promotes flexibility and code reuse.
 * 
 * 2. IMMUTABLE OBJECT PATTERN: Like Point3D, Line3D is immutable, ensuring thread safety and 
 *    predictable behavior in concurrent environments.
 * 
 * 3. FACTORY METHOD PATTERN: Multiple static factory methods provide clear, readable ways to create lines
 *    from different representations (points, vectors, parametric form).
 * 
 * 4. STRATEGY PATTERN: Different distance calculation methods implement various geometric algorithms
 *    (point-to-line, line-to-line, skew line distance).
 * 
 * 5. NULL OBJECT PATTERN: Uses Optional<> for methods that may not have valid results, avoiding null returns.
 * 
 * 6. TEMPLATE METHOD PATTERN: Common validation and logging patterns are consistently applied.
 * 
 * Data Structure Principles Demonstrated:
 * - Encapsulation: Line properties are derived from private Point3D endpoints
 * - Composition over Inheritance: Uses Point3D objects rather than extending them
 * - Immutability: Thread-safe design prevents state corruption
 * - Memory Efficiency: Reuses Point3D objects and caches computed properties
 * 
 * Algorithm Principles Demonstrated:
 * - Computational Geometry: Line-line intersection, closest point algorithms
 * - Vector Mathematics: Direction vectors, parametric line equations
 * - Linear Algebra: Projection algorithms, orthogonal vector calculations
 * - Numerical Analysis: Handling degenerate cases and floating-point precision
 * - Optimization: Efficient distance calculations using vector projections
 * 
 * Mathematical Foundation:
 * A line in 3D space can be represented parametrically as: P(t) = P₀ + t * d
 * where P₀ is a point on the line, d is the direction vector, and t is the parameter.
 * 
 * @author Spring Framework Style Documentation
 * @version 1.0
 */
public final class Line3D {
    
    private static final Logger logger = Logger.getLogger(Line3D.class.getName());
    
    // Tolerance for floating-point comparisons - demonstrates NUMERICAL STABILITY pattern
    private static final double EPSILON = 1e-10;
    
    // Immutable endpoints - demonstrates COMPOSITION and VALUE OBJECT patterns
    private final Point3D startPoint;
    private final Point3D endPoint;
    
    // Lazily computed cached properties - demonstrates LAZY INITIALIZATION pattern
    private transient Point3D directionVector;
    private transient Double length;
    
    /**
     * Constructs a new Line3D between two specified points.
     * 
     * This constructor demonstrates the DEFENSIVE PROGRAMMING pattern by validating inputs
     * and ensuring the line is not degenerate (zero length). It follows Spring Framework
     * conventions of explicit parameter validation and clear error messaging.
     * 
     * The constructor enforces the mathematical requirement that a line must have two
     * distinct points, preventing degenerate cases that would cause division by zero
     * in subsequent calculations.
     * 
     * @param startPoint the starting point of the line segment (must not be null)
     * @param endPoint the ending point of the line segment (must not be null and different from start)
     * @throws IllegalArgumentException if either point is null or points are identical
     */
    public Line3D(Point3D startPoint, Point3D endPoint) {
        if (startPoint == null) {
            logger.log(Level.SEVERE, "Cannot create Line3D: start point is null");
            throw new IllegalArgumentException("Start point cannot be null");
        }
        
        if (endPoint == null) {
            logger.log(Level.SEVERE, "Cannot create Line3D: end point is null");
            throw new IllegalArgumentException("End point cannot be null");
        }
        
        if (startPoint.equals(endPoint)) {
            logger.log(Level.SEVERE, "Cannot create Line3D: start and end points are identical: {0}", startPoint);
            throw new IllegalArgumentException("Start and end points must be different");
        }
        
        this.startPoint = startPoint;
        this.endPoint = endPoint;
        
        logger.log(Level.INFO, "Created new Line3D from {0} to {1}", new Object[]{startPoint, endPoint});
    }
    
    /**
     * Static factory method to create a Line3D from a point and direction vector.
     * 
     * This method demonstrates the FACTORY METHOD pattern, providing an intuitive way
     * to create lines using point-direction representation, which is common in
     * computational geometry and computer graphics.
     * 
     * The method converts from parametric form P(t) = P₀ + t*d to the two-point
     * representation used internally, using unit length for the direction vector.
     * 
     * This approach follows Spring's philosophy of providing multiple factory methods
     * with descriptive names rather than overloaded constructors.
     * 
     * @param point a point on the line (must not be null)
     * @param direction the direction vector (must not be null or zero vector)
     * @param length the length of the line segment from the point
     * @return a new Line3D representing the specified line segment
     * @throws IllegalArgumentException if parameters are invalid
     */
    public static Line3D fromPointAndDirection(Point3D point, Point3D direction, double length) {
        if (point == null) {
            logger.log(Level.SEVERE, "Cannot create Line3D: point is null");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        if (direction == null) {
            logger.log(Level.SEVERE, "Cannot create Line3D: direction vector is null");
            throw new IllegalArgumentException("Direction vector cannot be null");
        }
        
        if (length <= 0) {
            logger.log(Level.SEVERE, "Cannot create Line3D: invalid length {0}", length);
            throw new IllegalArgumentException("Length must be positive");
        }
        
        // Normalize direction and scale by length
        Point3D normalizedDirection = direction.normalize();
        Point3D endPoint = point.translate(
            normalizedDirection.getX() * length,
            normalizedDirection.getY() * length,
            normalizedDirection.getZ() * length
        );
        
        logger.log(Level.INFO, "Created Line3D from point {0} with direction {1} and length {2}", 
                  new Object[]{point, direction, length});
        
        return new Line3D(point, endPoint);
    }
    
    /**
     * Static factory method to create a Line3D representing an axis-aligned line.
     * 
     * This convenience method demonstrates the CONVENIENCE FACTORY pattern, providing
     * easy creation of common geometric primitives. Axis-aligned lines are frequently
     * used in 3D graphics, collision detection, and spatial indexing algorithms.
     * 
     * @param origin the starting point of the axis line
     * @param axis the axis direction ('X', 'Y', or 'Z')
     * @param length the length of the line segment
     * @return a new Line3D along the specified axis
     * @throws IllegalArgumentException if parameters are invalid
     */
    public static Line3D createAxisAligned(Point3D origin, char axis, double length) {
        if (origin == null) {
            logger.log(Level.SEVERE, "Cannot create axis-aligned line: origin is null");
            throw new IllegalArgumentException("Origin cannot be null");
        }
        
        if (length <= 0) {
            logger.log(Level.SEVERE, "Cannot create axis-aligned line: invalid length {0}", length);
            throw new IllegalArgumentException("Length must be positive");
        }
        
        Point3D endPoint;
        switch (Character.toLowerCase(axis)) {
            case 'x':
                endPoint = origin.translate(length, 0, 0);
                break;
            case 'y':
                endPoint = origin.translate(0, length, 0);
                break;
            case 'z':
                endPoint = origin.translate(0, 0, length);
                break;
            default:
                logger.log(Level.SEVERE, "Invalid axis specified: {0}", axis);
                throw new IllegalArgumentException("Axis must be 'X', 'Y', or 'Z'");
        }
        
        logger.log(Level.INFO, "Created axis-aligned Line3D along {0}-axis from {1} with length {2}", 
                  new Object[]{axis, origin, length});
        
        return new Line3D(origin, endPoint);
    }
    
    /**
     * Calculates the length of this line segment.
     * 
     * This method implements the LAZY INITIALIZATION pattern, computing and caching
     * the length on first access. Subsequent calls return the cached value, demonstrating
     * an optimization technique for expensive computations.
     * 
     * The length calculation uses the Euclidean distance formula, leveraging the
     * distanceTo() method from our Point3D class, which demonstrates the benefits
     * of the COMPOSITION pattern - reusing well-tested functionality.
     * 
     * Time complexity: O(1) after first computation due to caching
     * Space complexity: O(1) - single cached double value
     * 
     * @return the length of the line segment
     */
    public double length() {
        if (length == null) {
            length = startPoint.distanceTo(endPoint);
            logger.log(Level.INFO, "Calculated and cached length of line {0}: {1}", new Object[]{this, length});
        }
        return length;
    }
    
    /**
     * Returns the direction vector of this line.
     * 
     * The direction vector is computed as (endPoint - startPoint) and represents
     * the vector from start to end. This method also uses LAZY INITIALIZATION
     * for performance optimization.
     * 
     * The direction vector is fundamental for many geometric algorithms including
     * line intersection, closest point calculations, and parametric equations.
     * 
     * @return a Point3D representing the direction vector (not normalized)
     */
    public Point3D getDirectionVector() {
        if (directionVector == null) {
            directionVector = endPoint.translate(-startPoint.getX(), -startPoint.getY(), -startPoint.getZ());
            logger.log(Level.INFO, "Calculated and cached direction vector for line {0}: {1}", 
                      new Object[]{this, directionVector});
        }
        return directionVector;
    }
    
    /**
     * Returns the normalized (unit) direction vector of this line.
     * 
     * The unit direction vector has length 1 and points in the same direction as
     * the line. This is essential for many geometric calculations where direction
     * matters but magnitude doesn't.
     * 
     * @return a Point3D representing the normalized direction vector
     */
    public Point3D getUnitDirectionVector() {
        Point3D direction = getDirectionVector();
        Point3D unitDirection = direction.normalize();
        
        logger.log(Level.INFO, "Calculated unit direction vector for line {0}: {1}", 
                  new Object[]{this, unitDirection});
        
        return unitDirection;
    }
    
    /**
     * Calculates the midpoint of this line segment.
     * 
     * The midpoint is computed using the formula: midpoint = (start + end) / 2
     * This demonstrates basic vector arithmetic and is useful for line subdivision,
     * bounding box calculations, and geometric center computations.
     * 
     * @return a Point3D representing the midpoint of the line segment
     */
    public Point3D getMidpoint() {
        Point3D midpoint = new Point3D(
            (startPoint.getX() + endPoint.getX()) / 2.0,
            (startPoint.getY() + endPoint.getY()) / 2.0,
            (startPoint.getZ() + endPoint.getZ()) / 2.0
        );
        
        logger.log(Level.INFO, "Calculated midpoint of line {0}: {1}", new Object[]{this, midpoint});
        
        return midpoint;
    }
    
    /**
     * Calculates the shortest distance from a point to this line segment.
     * 
     * This method implements a fundamental computational geometry algorithm with three cases:
     * 1. Point projects onto the line segment: distance to projection point
     * 2. Point projects before start: distance to start point  
     * 3. Point projects after end: distance to end point
     * 
     * The algorithm uses vector projection: proj = (AP · AB) / |AB|²
     * where A is start point, B is end point, P is the query point.
     * 
     * This demonstrates the ALGORITHM STRATEGY pattern, selecting the appropriate
     * calculation method based on the geometric configuration.
     * 
     * Time complexity: O(1) - constant time vector operations
     * Space complexity: O(1) - temporary variables only
     * 
     * @param point the point to calculate distance from (must not be null)
     * @return the shortest distance from the point to this line segment
     * @throws IllegalArgumentException if point is null
     */
    public double distanceToPoint(Point3D point) {
        if (point == null) {
            logger.log(Level.SEVERE, "Cannot calculate distance: point is null");
            throw new IllegalArgumentException("Point cannot be null");
        }
        
        Point3D direction = getDirectionVector();
        Point3D startToPoint = point.translate(-startPoint.getX(), -startPoint.getY(), -startPoint.getZ());
        
        // Calculate projection parameter t
        double dot = startToPoint.dotProduct(direction);
        double lengthSquared = direction.magnitude() * direction.magnitude();
        double t = dot / lengthSquared;
        
        Point3D closestPoint;
        double distance;
        
        if (t < 0.0) {
            // Point projects before line start
            closestPoint = startPoint;
            distance = point.distanceTo(startPoint);
            logger.log(Level.INFO, "Point {0} projects before line start, closest point is start: {1}", 
                      new Object[]{point, startPoint});
        } else if (t > 1.0) {
            // Point projects after line end
            closestPoint = endPoint;
            distance = point.distanceTo(endPoint);
            logger.log(Level.INFO, "Point {0} projects after line end, closest point is end: {1}", 
                      new Object[]{point, endPoint});
        } else {
            // Point projects onto line segment
            closestPoint = new Point3D(
                startPoint.getX() + t * direction.getX(),
                startPoint.getY() + t * direction.getY(),
                startPoint.getZ() + t * direction.getZ()
            );
            distance = point.distanceTo(closestPoint);
            logger.log(Level.INFO, "Point {0} projects onto line segment, closest point: {1}", 
                      new Object[]{point, closestPoint});
        }
        
        logger.log(Level.INFO, "Distance from point {0} to line {1}: {2}", 
                  new Object[]{point, this, distance});
        
        return distance;
    }
    
    /**
     * Calculates the shortest distance between this line and another line in 3D space.
     * 
     * This method implements one of the most complex algorithms in computational geometry,
     * handling all possible cases for two lines in 3D space:
     * 
     * 1. PARALLEL LINES: Use point-to-line distance
     * 2. INTERSECTING LINES: Distance is zero at intersection point
     * 3. SKEW LINES: Use the formula involving cross and dot products
     * 
     * For skew lines, the distance formula is:
     * distance = |((P₁-P₂) · (d₁ × d₂))| / |d₁ × d₂|
     * 
     * where P₁, P₂ are points on the lines and d₁, d₂ are direction vectors.
     * 
     * This demonstrates the TEMPLATE METHOD pattern with different strategies
     * for different geometric configurations, and shows advanced vector mathematics.
     * 
     * Time complexity: O(1) - constant time vector operations
     * Space complexity: O(1) - temporary calculation variables
     * 
     * @param other the other line to calculate distance to (must not be null)
     * @return the shortest distance between the two line segments
     * @throws IllegalArgumentException if other line is null
     */
    public double shortestDistanceToLine(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot calculate distance: other line is null");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        Point3D d1 = this.getDirectionVector();
        Point3D d2 = other.getDirectionVector();
        Point3D crossProduct = d1.crossProduct(d2);
        
        // Check if lines are parallel (cross product magnitude ≈ 0)
        double crossMagnitude = crossProduct.magnitude();
        if (crossMagnitude < EPSILON) {
            logger.log(Level.INFO, "Lines are parallel, using point-to-line distance");
            return other.distanceToPoint(this.startPoint);
        }
        
        // Lines are not parallel - calculate skew line distance
        Point3D startDifference = other.startPoint.translate(
            -this.startPoint.getX(),
            -this.startPoint.getY(), 
            -this.startPoint.getZ()
        );
        
        double distance = Math.abs(startDifference.dotProduct(crossProduct)) / crossMagnitude;
        
        logger.log(Level.INFO, "Calculated shortest distance between skew lines {0} and {1}: {2}", 
                  new Object[]{this, other, distance});
        
        return distance;
    }
    
    /**
     * Finds the intersection point of this line with another line, if it exists.
     * 
     * This method solves the system of parametric equations:
     * P₁(s) = A₁ + s * d₁
     * P₂(t) = A₂ + t * d₂
     * 
     * Setting P₁(s) = P₂(t) and solving for s and t. If a solution exists within
     * the parameter ranges [0,1] for both lines, the lines intersect.
     * 
     * This demonstrates the NULL OBJECT pattern using Optional<> to handle the case
     * where no intersection exists, avoiding null returns and potential NullPointerExceptions.
     * 
     * The method implements sophisticated numerical analysis to handle floating-point
     * precision issues and degenerate cases.
     * 
     * @param other the other line to find intersection with
     * @return Optional containing intersection point if it exists, empty otherwise
     * @throws IllegalArgumentException if other line is null
     */
    public Optional<Point3D> findIntersection(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot find intersection: other line is null");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        Point3D d1 = this.getDirectionVector();
        Point3D d2 = other.getDirectionVector();
        Point3D startDiff = other.startPoint.translate(
            -this.startPoint.getX(),
            -this.startPoint.getY(),
            -this.startPoint.getZ()
        );
        
        // Check if lines are parallel
        Point3D crossProduct = d1.crossProduct(d2);
        if (crossProduct.magnitude() < EPSILON) {
            logger.log(Level.WARNING, "Lines are parallel, no unique intersection point");
            return Optional.empty();
        }
        
        // Solve the system of equations using cross products and dot products
        // This is a simplified approach - in practice, you'd use more robust linear algebra
        Point3D cross1 = startDiff.crossProduct(d2);
        Point3D cross2 = startDiff.crossProduct(d1);
        
        double t1 = cross1.dotProduct(crossProduct) / crossProduct.dotProduct(crossProduct);
        double t2 = cross2.dotProduct(crossProduct) / crossProduct.dotProduct(crossProduct);
        
        // Check if intersection occurs within both line segments
        if (t1 >= 0 && t1 <= 1 && t2 >= 0 && t2 <= 1) {
            Point3D intersection = new Point3D(
                this.startPoint.getX() + t1 * d1.getX(),
                this.startPoint.getY() + t1 * d1.getY(),
                this.startPoint.getZ() + t1 * d1.getZ()
            );
            
            logger.log(Level.INFO, "Found intersection between lines {0} and {1} at point {2}", 
                      new Object[]{this, other, intersection});
            
            return Optional.of(intersection);
        }
        
        logger.log(Level.INFO, "Lines {0} and {1} do not intersect within their segments", 
                  new Object[]{this, other});
        
        return Optional.empty();
    }
    
    /**
     * Determines if this line is parallel to another line.
     * 
     * Two lines are parallel if their direction vectors are parallel (cross product ≈ 0).
     * This method uses a tolerance-based comparison to handle floating-point precision issues.
     * 
     * This demonstrates the TOLERANCE-BASED EQUALITY pattern for floating-point comparisons,
     * which is essential in computational geometry where exact equality rarely occurs.
     * 
     * @param other the other line to check parallelism with
     * @return true if the lines are parallel, false otherwise
     * @throws IllegalArgumentException if other line is null
     */
    public boolean isParallelTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot check parallelism: other line is null");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        Point3D d1 = this.getUnitDirectionVector();
        Point3D d2 = other.getUnitDirectionVector();
        Point3D crossProduct = d1.crossProduct(d2);
        
        boolean parallel = crossProduct.magnitude() < EPSILON;
        
        logger.log(Level.INFO, "Lines {0} and {1} are {2}parallel", 
                  new Object[]{this, other, parallel ? "" : "not "});
        
        return parallel;
    }
    
    /**
     * Determines if this line is perpendicular to another line.
     * 
     * Two lines are perpendicular if their direction vectors are orthogonal (dot product ≈ 0).
     * This is another example of tolerance-based floating-point comparison.
     * 
     * @param other the other line to check perpendicularity with
     * @return true if the lines are perpendicular, false otherwise
     * @throws IllegalArgumentException if other line is null
     */
    public boolean isPerpendicularTo(Line3D other) {
        if (other == null) {
            logger.log(Level.SEVERE, "Cannot check perpendicularity: other line is null");
            throw new IllegalArgumentException("Other line cannot be null");
        }
        
        Point3D d1 = this.getUnitDirectionVector();
        Point3D d2 = other.getUnitDirectionVector();
        double dotProduct = Math.abs(d1.dotProduct(d2));
        
        boolean perpendicular = dotProduct < EPSILON;
        
        logger.log(Level.INFO, "Lines {0} and {1} are {2}perpendicular", 
                  new Object[]{this, other, perpendicular ? "" : "not "});
        
        return perpendicular;
    }
    
    /**
     * Returns a point on this line at the specified parameter value.
     * 
     * This method implements the parametric line equation: P(t) = P₀ + t * d
     * where t=0 gives the start point, t=1 gives the end point, and values
     * outside [0,1] extend beyond the line segment.
     * 
     * This is useful for line subdivision, animation along paths, and
     * implementing various geometric algorithms.
     * 
     * @param parameter the parameter value (0.0 = start point, 1.0 = end point)
     * @return a Point3D on the line at the specified parameter
     */
    public Point3D getPointAt(double parameter) {
        if (Double.isNaN(parameter) || Double.isInfinite(parameter)) {
            logger.log(Level.WARNING, "Invalid parameter value: {0}, returning start point", parameter);
            return startPoint;
        }
        
        Point3D direction = getDirectionVector();
        Point3D result = new Point3D(
            startPoint.getX() + parameter * direction.getX(),
            startPoint.getY() + parameter * direction.getY(),
            startPoint.getZ() + parameter * direction.getZ()
        );
        
        logger.log(Level.INFO, "Calculated point at parameter {0} on line {1}: {2}", 
                  new Object[]{parameter, this, result});
        
        return result;
    }
    
    // Accessor methods - demonstrate ENCAPSULATION principle
    
    /**
     * Returns the start point of this line segment.
     * @return the start point
     */
    public Point3D getStartPoint() { return startPoint; }
    
    /**
     * Returns the end point of this line segment.  
     * @return the end point
     */
    public Point3D getEndPoint() { return endPoint; }
    
    /**
     * Indicates whether some other object is "equal to" this one.
     * 
     * Two Line3D objects are equal if they have the same start and end points.
     * This implementation handles the case where lines have opposite directions
     * but represent the same geometric line segment.
     * 
     * @param obj the reference object with which to compare
     * @return true if this object equals the obj argument; false otherwise
     */
    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        
        Line3D line3D = (Line3D) obj;
        
        // Check both orientations (A->B equals B->A for line segments)
        return (startPoint.equals(line3D.startPoint) && endPoint.equals(line3D.endPoint)) ||
               (startPoint.equals(line3D.endPoint) && endPoint.equals(line3D.startPoint));
    }
    
    /**
     * Returns a hash code value for this line.
     * 
     * The hash code is computed using both endpoints in a way that's orientation-independent,
     * ensuring that lines with swapped endpoints have the same hash code.
     * 
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        // Use XOR to make hash code orientation-independent
        return startPoint.hashCode() ^ endPoint.hashCode();
    }
    
    /**
     * Returns a string representation of this line.
     * 
     * The format clearly shows the start and end points with an arrow indicating direction.
     * This aids in debugging and logging geometric operations.
     * 
     * @return a string representation of the line
     */
    @Override
    public String toString() {
        return String.format("Line3D[%s -> %s]", startPoint, endPoint);
    }
}
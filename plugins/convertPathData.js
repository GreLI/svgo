'use strict';

exports.type = 'perItem';

exports.active = true;

exports.params = {
    applyTransforms: true,
    applyTransformsStroked: true,
    straightCurves: true,
    lineShorthands: true,
    curveSmoothShorthands: true,
    floatPrecision: 3,
    removeUseless: true,
    collapseRepeated: true,
    leadingZero: true,
    negativeExtraSpace: true,
    removeSmallCurves: false // 1 is a good value to start
};

var pathElems = require('./_collections.js').pathElems,
    path2js = require('./_path.js').path2js,
    js2path = require('./_path.js').js2path,
    applyTransforms = require('./_path.js').applyTransforms,
    hasMarkerMid;

/**
 * Convert absolute Path to relative,
 * collapse repeated instructions,
 * detect and convert Lineto shorthands,
 * remove useless instructions like "l0,0",
 * trim useless delimiters and leading zeros,
 * decrease accuracy of floating-point numbers.
 *
 * @see http://www.w3.org/TR/SVG/paths.html#PathData
 *
 * @param {Object} item current iteration item
 * @param {Object} params plugin params
 * @return {Boolean} if false, item will be filtered out
 *
 * @author Kir Belevich
 */
exports.fn = function(item, params) {

    if (item.isElem(pathElems) && item.hasAttr('d')) {

        hasMarkerMid = item.hasAttr('marker-mid');

        var data = path2js(item.attr('d').value);


        // TODO: get rid of functions returns
        if (data.length) {
            data = convertToRelative(data);

            if (params.applyTransforms) {
                data = applyTransforms(item, data, params.applyTransformsStroked, params.floatPrecision);
            }

            data = filters(data, params);

            if (params.collapseRepeated) {
                data = collapseRepeated(data, params);
            }

            item.pathJS = data;

            item.attr('d').value = js2path(data, params);
        }

    }

};

/**
 * Convert absolute path data coordinates to relative.
 *
 * @param {Array} path input path data
 * @param {Object} params plugin params
 * @return {Array} output path data
 */
function convertToRelative(path) {

    var instruction,
        data,
        newPoint,
        point = [0, 0],
        subpathPoint = [0, 0],
        index = 0,
        mM = false;

    path.forEach(function(item) {

        instruction = item.instruction;
        data = item.data;

        index++;

        // data !== !z
        if (data) {

            // already relative
            // recalculate current point
            if ('mcslqta'.indexOf(instruction) > -1) {

                newPoint = data.slice(-2);

                point[0] += newPoint[0];
                point[1] += newPoint[1];

                if (instruction === 'm') {
                    if (index === 1) {
                        instruction = 'M';
                        mM = true;
                    }

                    subpathPoint = point.slice(-2);
                }

            } else if (instruction === 'h') {

                point[0] += data[0];

            } else if (instruction === 'v') {

                point[1] += data[0];

            }

            // convert absolute path data coordinates to relative
            // M → m
            if (instruction === 'M') {

                if (index > 1) {
                    instruction = 'm';
                }

                // if "M" was not transformed from "m"
                if (!mM) {
                    data[0] -= point[0];
                    data[1] -= point[1];

                    point[0] += data[0];
                    point[1] += data[1];

                    subpathPoint = point.slice(-2);
                }

            }

            // L → l
            // T → t
            else if ('LT'.indexOf(instruction) > -1) {

                instruction = instruction.toLowerCase();

                // x y
                // 0 1
                data[0] -= point[0];
                data[1] -= point[1];

                point[0] += data[0];
                point[1] += data[1];

            // C → c
            } else if (instruction === 'C') {

                instruction = 'c';

                // x1 y1 x2 y2 x y
                // 0  1  2  3  4 5
                data[0] -= point[0];
                data[1] -= point[1];
                data[2] -= point[0];
                data[3] -= point[1];
                data[4] -= point[0];
                data[5] -= point[1];

                point[0] += data[4];
                point[1] += data[5];

            // S → s
            // Q → q
            } else if ('SQ'.indexOf(instruction) > -1) {

                instruction = instruction.toLowerCase();

                // x1 y1 x y
                // 0  1  2 3
                data[0] -= point[0];
                data[1] -= point[1];
                data[2] -= point[0];
                data[3] -= point[1];

                point[0] += data[2];
                point[1] += data[3];

            // A → a
            } else if (instruction === 'A') {

                instruction = 'a';

                // rx ry x-axis-rotation large-arc-flag sweep-flag x y
                // 0  1  2               3              4          5 6
                data[5] -= point[0];
                data[6] -= point[1];

                point[0] += data[5];
                point[1] += data[6];

            // H → h
            } else if (instruction === 'H') {

                instruction = 'h';

                data[0] -= point[0];

                point[0] += data[0];

            // V → v
            } else if (instruction === 'V') {

                instruction = 'v';

                data[0] -= point[1];

                point[1] += data[0];

            }

            item.instruction = instruction;
            item.data = data;
            item.absCoords = point.slice(-2);

        }

        // !data === z, reset current point
        else {
            point = subpathPoint;
            mM = false;
        }

    });

    return path;

}

/**
 * Main filters loop.
 *
 * @param {Array} path input path data
 * @param {Object} params plugin params
 * @return {Array} output path data
 */
function filters(path, params) {

    var instruction,
        data,
        prev,
        subpathPoint;

    path = path.filter(function(item, index, path) {

        instruction = item.instruction;
        data = item.data;

        if ('mM'.indexOf(instruction) > -1) subpathPoint = item.absCoords;

        if (data) {

            if (params.floatPrecision !== false) {
                data = roundData(data, params.floatPrecision);
            }

            // convert straight curves into lines segments
            if (params.straightCurves) {

                // c
                if (
                    instruction === 'c' &&
                    isCurveStraightLine(
                        [ 0, data[0], data[2], data[4] ],
                        [ 0, data[1], data[3], data[5] ]
                    )
                ) {
                    instruction = 'l';
                    data = data.slice(-2);
                }

                // s
                else if (instruction === 's') {

                    if (
                        isCurveStraightLine(
                            [ 0, data[0], data[2] ],
                            [ 0, data[1], data[3] ]
                        )
                    ) {
                        instruction = 'l';
                        data = data.slice(-2);
                    }

                }

                // q
                else if (
                    prev &&
                    instruction === 'q' &&
                    isCurveStraightLine(
                        [ 0, data[0], data[2] ],
                        [ 0, data[1], data[3] ]
                    )
                ) {
                    // save the original one for the future potential q + t conversion
                    item.original = {
                        instruction: instruction,
                        data: data
                    };

                    instruction = 'l';
                    data = data.slice(-2);
                }

                else if (instruction === 't') {

                    // q (original) + t
                    if (
                        prev &&
                        prev.original &&
                        prev.original.instruction === 'q'
                    ) {
                        if (isCurveStraightLine(
                            [ prev.original.data[0], prev.original.data[2], data[0] ],
                            [ prev.original.data[1], prev.original.data[3], data[1] ]
                        )) {
                            instruction = 'l';
                            data = data.slice(-2);
                        } else {
                            prev.instruction = 'q';
                            prev.data = prev.original.data;
                        }
                    }

                    // [^qt] + t
                    else if (!prev || 'qt'.indexOf(prev.instruction) === -1) {
                        instruction = 'l';
                        data = data.slice(-2);
                    }

                }

                // a
                else if (
                    instruction === 'a' &&
                    (data[0] === 0 || data[1] === 0)
                ) {
                    instruction = 'l';
                    data = data.slice(-2);
                }
            }

            // horizontal and vertical line shorthands
            // l 50 0 → h 50
            // l 0 50 → v 50
            if (
                params.lineShorthands &&
                instruction === 'l'
            ) {
                if (data[1] === 0) {
                    instruction = 'h';
                    data = [data[0]];
                } else if (data[0] === 0) {
                    instruction = 'v';
                    data = [data[1]];
                }
            }

            // convert curves into smooth shorthands
            if (params.curveSmoothShorthands && prev) {

                // curveto
                if (instruction === 'c') {

                    // c + c → c + s
                    if (
                        prev.instruction === 'c' &&
                        data[0] === -(prev.data[2] - prev.data[4]) &&
                        data[1] === -(prev.data[3] - prev.data[5])
                    ) {
                        instruction = 's';
                        data = data.slice(2);
                    }

                    // s + c → s + s
                    else if (
                        prev.instruction === 's' &&
                        data[0] === -(prev.data[0] - prev.data[2]) &&
                        data[1] === -(prev.data[1] - prev.data[3])
                    ) {
                        instruction = 's';
                        data = data.slice(2);
                    }

                    // [^cs] + c → [^cs] + s
                    else if (
                        'cs'.indexOf(prev.instruction) === -1 &&
                        data[0] === 0 &&
                        data[1] === 0
                    ) {
                        instruction = 's';
                        data = data.slice(2);
                    }

                }

                // quadratic Bézier curveto
                else if (instruction === 'q') {

                    // q + q → q + t
                    if (
                        prev.instruction === 'q' &&
                        data[0] === (prev.data[2] - prev.data[0]) &&
                        data[1] === (prev.data[3] - prev.data[1])
                    ) {
                        instruction = 't';
                        data = data.slice(2);
                    }

                    // t + q → t + t
                    else if (
                        prev.instruction === 't' &&
                        data[2] === prev.data[0] &&
                        data[3] === prev.data[1]
                    ) {
                        instruction = 't';
                        data = data.slice(2);
                    }

                }

            }

            // remove useless non-first path segments
            if (params.removeUseless) {

                // m 0,0 / l 0,0 / h 0 / v 0 / q 0,0 0,0 / t 0,0 / c 0,0 0,0 0,0 / s 0,0 0,0
                if (
                    (
                     'lhvqtcs'.indexOf(instruction) > -1
                    ) &&
                    data.every(function(i) { return i === 0; })
                ) {
                    return false;
                }

                // a 25,25 -30 0,1 0,0
                if (
                    instruction === 'a' &&
                    data[5] === 0 &&
                    data[6] === 0
                ) {
                    return false;
                }

            }

            if (
                params.removeSmallCurves &&
                'qc'.indexOf(instruction) > -1 &&
                Math.pow(item.absCoords[0] - path[index - 1].absCoords[0], 2) +
                Math.pow(item.absCoords[1] - path[index - 1].absCoords[1], 2) <
                Math.pow(params.removeSmallCurves, 2)
            ) {
                var adjancentItems = getAdjancentItems(path, index),
                    controlPoints = [
                        data[0],
                        data[1],
                        instruction === 'c' ? data[2] - data[4] : data[0] - data[2],
                        instruction === 'c' ? data[3] - data[5] : data[1] - data[3],
                    ];

                // Check if the curve is relatively blunt (angle > 45°).
                if (adjancentItems && checkAngles(controlPoints)) {
                    var items = adjancentItems.items,
                        linesPoints = adjancentItems.points,
                        nextM = adjancentItems.nextM;
                    // Get lines intersection point. The angle check ensures they aren't parallel.
                    var crossPoint = getCrossingPoint(linesPoints),
                        newPoints = [
                            crossPoint[0] - linesPoints[0],
                            crossPoint[1] - linesPoints[1],
                            linesPoints[6] - crossPoint[0],
                            linesPoints[7] - crossPoint[1]
                        ];

                    if (params.floatPrecision) {
                        newPoints = roundData(newPoints, params.floatPrecision);
                    }

                    // Modify adjancent lines
                    for (var i = 0; i < items.length; i++) {
                        switch (items[i].instruction) {
                            case 'c':
                                items[i].data[4] = newPoints[2*i];
                                items[i].data[5] = newPoints[2*i + 1];
                            break;
                            case 'm':
                            case 'l': items[i].data[1] = newPoints[2*i + 1]; // falls through for x-coord
                            case 'h': items[i].data[0] = newPoints[2*i];
                            break;
                            case 'v': items[i].data[0] = newPoints[2*i + 1];
                        }
                    }
                    // Update prev and current item too, since it will stay along the filter process.
                    items[0].absCoords[0] = item.absCoords[0] = crossPoint[0];
                    items[0].absCoords[1] = item.absCoords[1] = crossPoint[1];

                    if (nextM && nextM.instruction === 'm') {
                        nextM.data[0] = nextM.absCoords[0] - subpathPoint[0];
                        nextM.data[1] = nextM.absCoords[1] - subpathPoint[1];
                    }

                    // Remove the curve
                    return false;
                }
            }

            item.instruction = instruction;
            item.data = data;

            prev = item;

        }

        return true;

    });

    return path;

}

/**
 * Collapse repeated instructions data
 *
 * @param {Array} path input path data
 * @return {Array} output path data
 */
function collapseRepeated(path) {

    var prev;

    path = path.filter(function(item) {

        if (
            !hasMarkerMid &&
            prev &&
            item.instruction === prev.instruction
        ) {
            // increase previous h or v data with current
            if ((item.instruction === 'h' || item.instruction === 'v') && (prev.data[0] >= 0) == (item.data[0] >= 0)) {
                prev.data[0] += item.data[0];
            // concat previous data with current if it is not z
            } else if (item.instruction !== 'z') {
                prev.data = prev.data.concat(item.data);
            }

            // filter current item
            return false;
        }

        prev = item;

        return true;

    });

    return path;

}

/**
 * Decrease accuracy of floating-point numbers
 * in path data keeping a specified number of decimals.
 *
 * @param {Array} data input data array
 * @param {Number} fixed number of decimals
 * @return {Array} output data array
 */
function roundData(data, fixed) {

    return data.map(function(num) {
        return +num.toFixed(fixed);
    });

}

/**
 * Checks if curve is a straight line by calculating a polygon area.
 *
 * @see http://www.mathopenref.com/coordpolygonarea2.html
 *
 * @param {Array} xs array of curve points x-coordinates
 * @param {Array} ys array of curve points y-coordinates
 * @return {Boolean}
 */

function isCurveStraightLine(xs, ys) {

    var points = xs.length,
        area = 0,
        j = points - 1;

    for (var i=0; i < points; i++) {
        area += (xs[j] + xs[i]) * (ys[j] - ys[i]);
        j = i;
    }

    if (+area.toFixed(2)) return false;

    return true;

}

/**
 * Check if angles more then 45° 
 *
 * @param {Array} coordinates of adjancent lines points (8 = 4×2)
 * @return {Array} coordinate of lines intersection point, false if none (parallel lines)
 */

function checkAngles(points) {

    if (
        !points.length ||
        points.length !== 4 ||
        !points.every(function(n){ return !isNaN(n) && isFinite(n) })
    )
        return false;

    var num = points[0]*points[1] + points[2]*points[3],
        denom1 = points[0]*points[0] + points[1]*points[1],
        denom2 = points[2]*points[2] + points[3]*points[3];

    // Check the angle between the curve tangents by computing cos^2.
    // (Square to avoid taking root but sign is preserved.)
    return (
        denom1 > 0 &&
        denom2 > 0 &&
        num * Math.abs(num) / (denom1*denom2) < .5
    ) ? true : false;

}

/**
 * Computes line equation from points and gets the intersection point.
 *
 * @param {Array} coordinates of adjancent lines points (8 numbers for 4 points)
 * @return {Array} coordinate of lines intersection point, false if none (parallel lines or bad arguments)
 */

function getCrossingPoint(points) {
    if (
        !points.length ||
        points.length !== 8 ||
        !points.every(function(n){ return !isNaN(n) && isFinite(n) })
    )
        return false;

    var // Prev line equation parameters.
        a = points[1] - points[3], // y1 - y2
        b = points[2] - points[0], // x2 - x1
        c = points[0]*points[3] - points[2]*points[1], // x1*y2 - x2*y1
        // Next line equation parameters
        a1 = points[5] - points[7], // y1 - y2
        b1 = points[6] - points[4], // x2 - x1
        c1 = points[4]*points[7] - points[5]*points[6], // x1*y2 - x2*y1
        denom1 = (a*b1 - a1*b),
        denom2 = (b*a1 - b1*a);

    return denom1 !== 0 && denom2 !== 0 ? [
        (c1*b - c*b1) / denom1,
        (c1*a - c*a1) / denom2
    ] : false;

}

function getAdjancentItems(path, index) {

    // Simplest case: prev and next items are lines.
    if (
        'lhv'.indexOf(path[index - 1].instruction) > -1 &&
        (
            'lhv'.indexOf(path[index + 1].instruction) > -1 ||
            path[index + 1].instruction === 'c' &&
            isCurveStraightLine(
                [ 0, path[index + 1].data[0], path[index + 1].data[2], path[index + 1].data[4] ],
                [ 0, path[index + 1].data[1], path[index + 1].data[3], path[index + 1].data[5] ]
            )
        )
    ) {
        return {
            items: [
                path[index - 1],
                path[index + 1]
            ],
            points: [].concat(
                path[index - 2].absCoords,
                path[index - 1].absCoords,
                path[index].absCoords,
                path[index + 1].absCoords
            )
        }
    }

    var eps = .01;

    // The curve at the start.
    if (
        'Mm'.indexOf(path[index - 1].instruction) > -1 &&
        (
            'lhv'.indexOf(path[index + 1].instruction) > -1 ||
            path[index + 1].instruction === 'c' &&
            isCurveStraightLine(
                [ 0, path[index + 1].data[0], path[index + 1].data[2], path[index + 1].data[4] ],
                [ 0, path[index + 1].data[1], path[index + 1].data[3], path[index + 1].data[5] ]
            )
        )
    ) {

        var lastItem,
            startPoint;

        // Find the last item.
        for (var i = index; path[++i] && !'mz'.indexOf(path[i].instruction); )
            lastItem = path[i];

        // Path is closed by additional line formed automatically.
        if (
            Math.abs(lastItem.absCoords[0] - path[index - 1][0]) > eps ||
            Math.abs(lastItem.absCoords[1] - path[index - 1][1]) > eps
        ) {
            startPoint = path[i - 1].absCoords
        }

        // Last item is a line and closes the path.
        else if (
            'lhv'.indexOf(lastItem.instruction) > -1 &&
            Math.abs(lastItem.absCoords[0] - path[index - 1][0]) < eps &&
            Math.abs(lastItem.absCoords[1] - path[index - 1][1]) < eps
        ) {
            startPoint = path[i - 2].absCoords;
        }
        if (!startPoint) return false;

        var adjancentItems = {
                items: [
                    path[index - 1],
                    path[index + 1]
                ],
                points: [].concat(
                    startPoint,
                    path[index - 1].absCoords,
                    path[index].absCoords,
                    path[index + 1].absCoords
                )
            },
            nextM;

        // Find next move command
        for (var j = i; (nextM = path[j]) && 'Mm'.indexOf(nextM.instruction) === -1; j++);
        if (nextM) adjancentItems.nextM = nextM;

        return adjancentItems;
    }

    // The curve at the end.
    if (
        'lhv'.indexOf(path[index - 1].instruction) > -1 &&
        (
            !path[index + 1] ||
            'mz'.indexOf(path[index + 1].instruction) > -1
        )
    ) {

        // Find the start.
        for (var i = index; path[--i] && 'Mm'.indexOf(path[i].instruction) === -1; );

        if (!path[i]) return false;

        var items = [ path[index - 1] ], // No need in updating next sibling.
            endPoint;

        // Path is closed by additional line formed automatically.
        if (
            Math.abs(path[index][1] - path[i].absCoords[0]) > eps ||
            Math.abs(path[index][0] - path[i].absCoords[1]) > eps
        ) {
            endPoint = path[i].absCoords;
        }

        // First item is a line.
        else if ('lhv'.indexOf(path[i + 1].instruction) > -1) {
            endPoint = path[i + 1].absCoords;
        }

        if (!endPoint)  return false;

        var adjancentItems = {
                items: items,
                points: [].concat(
                    path[index - 2].absCoords,
                    path[index - 1].absCoords,
                    path[index].absCoords,
                    endPoint
                )
            },
            nextM;

        // Find next move command
        for (var i = index + 1; (nextM = path[i]) && 'Mm'.indexOf(nextM.instruction) === -1; i++);
        if (nextM) adjancentItems.nextM = nextM;

        return adjancentItems;

    }

    return false;

}

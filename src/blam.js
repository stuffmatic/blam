/*
 *TODO:
 * - pass image URL and control point positions etc with url:
 *      lolz.com/?img=rofl.com&blblaimg.png&npts=2&horiz=0&vp1x=0.33&vp21y=0.24&...&horizx=0.5&horizy=0.4
 * - save params as a cookie
 * - limit float precision in query string
 

 
 */

//ids for dynamic text elements
var ID_CALIBRATION_INFO_LABEL = "calibrationInfoLabel";
var ID_INFO_LABEL = "infoLabel";
var ID_FOCAL_LENGTH_LABEL = "focalLengthLabel";
var ID_FIELD_OF_VIEW_LABEL = "fieldOfViewLabel";
var ID_AXIS_ANGLE_LABEL = "axisAngleLabel";
var ID_EULER_LABEL = "eulerLabel";
var ID_QUAT_LABEL = "quatLabel";
var ID_TRANSLATION_LABEL = "translationLabel";
var ID_MATRIX_LABEL = "matrixLabel";
var ID_SHARE_URL_LABEL = "shareURLLabel";

//ids for other UI elements
var ID_CAMERA_MENU = "cameraMenu";
var ID_IMAGE_CANVAS = "imageCanvas";
var ID_OVERLAY_CANVAS = "overlayCanvas";
var ID_SPLASH_CONTAINER = "splashContainer";
var ID_SENSOR_WIDTH_TEXTFIELD = "sensorWidth";
var ID_SENSOR_HEIGHT_TEXTFIELD = "sensorHeight";
var ID_HORIZON_CHECKBOX = "horizonCheckBox";
var ID_NUM_VANISHING_POINTS = "numVanishingPoints";

//colors
var COL_VLX = "#aa2200";
var COL_VLY = "#22aa00";
var COL_VLZ = "#2175aa";
var COL_HORIZ = "#aaaaaa";
var COL_PRINCIPAL = "#ffaa00";
var COL_ORIGIN = "#f0f0f0";

//parameter names
var PARAM_NUM_VP = "nvp";
var PARAM_HORIZ_START_X = "hstartx";
var PARAM_HORIZ_START_Y = "hstarty";
var PARAM_HORIZ_END_X = "hendx";
var PARAM_HORIZ_END_Y = "hendy";
var PARAM_HORIZ_ENABLED = "henabled";
var PARAM_CAMERA_PRESET = "preset";
var PARAM_SENSOR_WIDTH = "sensorw";
var PARAM_SENSOR_HEIGHT = "sensorh";

var PARAM_IMAGE_URL = "imageurl";

var PARAM_VLX1_START_X = "vlx1startx";
var PARAM_VLX1_START_Y = "vlx1starty";
var PARAM_VLX1_END_X = "vlx1endx";
var PARAM_VLX1_END_Y = "vlx1endy";

var PARAM_VLX2_START_X = "vlx2startx";
var PARAM_VLX2_START_Y = "vlx2starty";
var PARAM_VLX2_END_X = "vlx2endx";
var PARAM_VLX2_END_Y = "vlx2endy";

var PARAM_VLY1_START_X = "vly1startx";
var PARAM_VLY1_START_Y = "vly1starty";
var PARAM_VLY1_END_X = "vly1endx";
var PARAM_VLY1_END_Y = "vly1endy";

var PARAM_VLY2_START_X = "vly2startx";
var PARAM_VLY2_START_Y = "vly2starty";
var PARAM_VLY2_END_X = "vly2endx";
var PARAM_VLY2_END_Y = "vly2endy";

var PARAM_VLZ1_START_X = "vlz1startx";
var PARAM_VLZ1_START_Y = "vlz1starty";
var PARAM_VLZ1_END_X = "vlz1endx";
var PARAM_VLZ1_END_Y = "vlz1endy";

var PARAM_VLZ2_START_X = "vlz2startx";
var PARAM_VLZ2_START_Y = "vlz2starty";
var PARAM_VLZ2_END_X = "vlz2endx";
var PARAM_VLZ2_END_Y = "vlz2endy";

var PARAM_ORIGIN_X = "originx";
var PARAM_ORIGIN_Y = "originy";

var PARAM_PRINCIPAL_X = "principalx";
var PARAM_PRINCIPAL_Y = "principaly";

var DEFAULT_STATE_QUERY_STRING =
"nvp=1&hstartx=0.07064765348896995&hstarty=0.6594830707960123&hendx=0.8727585863977095&hendy=0.6592779587049186&preset=0&sensorw=36&sensorh=24&vlx1startx=0.7233477088948785&vlx1starty=0.34374722970949384&vlx1endx=0.4334775766642072&vlx1endy=0.036275652596892805&vlx2startx=0.23675727210338274&vlx2starty=0.7317300285559857&vlx2endx=0.5164880454054931&vlx2endy=0.8748507371700606&vly1startx=0.49564463906230743&vly1starty=0.04225672585462735&vly1endx=0.25520158223251077&vly1endy=0.32984649130726973&vly2startx=0.7352133883758031&vly2starty=0.7245420231157701&vly2endx=0.44814967560232144&vly2endy=0.8723733049961717&vlz1startx=0.25957836207756724&vlz1starty=0.7912823007673254&vlz1endx=0.2918740340146725&vlz1endy=0.24388163777628594&vlz2startx=0.6695117999548136&vlz2starty=0.22641578288706898&vlz2endx=0.6921403652570008&vlz2endy=0.8033319447839186&originx=0.459505893853026&originy=0.5962724073192436&imageurl=webcube.png&henabled=0&principalx=0.5028573128555447&principaly=0.521936657940252";
//a dictionary defining the state of the app.
//used for loading and saving.
var appState =
{
    
};


/*
 *control point (CP) indices:
 */
//vanishing lines
var CP_VLX1_START = 0;
var CP_VLX1_END = 1;

var CP_VLX2_START = 2;
var CP_VLX2_END = 3;

var CP_VLY1_START = 4;
var CP_VLY1_END = 5;

var CP_VLY2_START = 6;
var CP_VLY2_END = 7;

var CP_VLZ1_START = 8;
var CP_VLZ1_END = 9;

var CP_VLZ2_START = 10;
var CP_VLZ2_END = 11;

//horizon line
var CP_HORIZ_START = 12;
var CP_HORIZ_END = 13;
//principal point
var CP_PRINCIPAL = 14;
//origin
var CP_ORIGIN = 15;
//total number of points
var CP_COUNT = 16;

//control points (vanishing lines, horizon, principal point, etc)
var controlPoints =
[
    //x
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    //y
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    //z
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    //horizon line
    [Math.random(), Math.random()],
    [Math.random(), Math.random()],
    
    //principal point
    [Math.random(), Math.random()],
    
    //origin
    [Math.random(), Math.random()]
];

//
var vanishingPoints = [[0, 0], [0, 0], [0, 0]];

//various constants
var SELECTION_RADIUS = 8;

/**
 *
 */
function length(vec)
{
    var lSq = 0.0;
    
    for (var i = 0; i < vec.length; i++)
    {
        lSq += vec[i] * vec[i];
    }
    
    return Math.sqrt(lSq);    
}

/**
 *
 */
function normalize(vec)
{
    var l = length(vec);
    
    if (l == 0.0)
    {
        return;
    }
    
    for (var i = 0; i < vec.length; i++)
    {
        vec[i] /= l;
    }
    
    return vec;
}

/**
 *
 */
function dot(v1, v2)
{
    var ret = 0.0;
    
    for (var i = 0; i < v1.length; i++)
    {
        ret += v1[i] * v2[i];
    }
    
    return ret;
}

function vecSubtract(v1, v2)
{
    console.log("subtr: v1 = " + v1 + ", v2 = " + v2);
    return [v2[0] - v1[0], v2[1] - v1[1]];
}

function truncateFloat(f)
{
    var i = parseInt(f * 1000);
    return i * 0.001;
}

/**
 * Computes the coordinates of the second vanishing point
 * based on the first, a focal length, the center of projection and 
 * the desired horizon tilt angle. The equations here are derived from 
 * section 3.2 "Determining the focal length from a single image".
 * @param Fu the first vanishing point in normalized image coordinates.
 * @param f the relative focal length.
 * @param P the center of projection in normalized image coordinates.
 * @param horizonDir The desired horizon direction.
 * @return The coordinates of the second vanishing point.
 * @see http://www.irisa.fr/prive/kadi/Reconstruction/paper.ps.gz
 */
function computeSecondVanishingPoint(Fu, f, P, horizonDir)
{
    //find the second vanishing point
    //TODO 1: take principal point into account here
    //TODO 2: if the first vanishing point coincides with the image center,
    //        these lines won't work, but this case should be handled somehow.
    var k = -(Fu[0] * Fu[0] + Fu[1] * Fu[1] + f * f) / (Fu[0] * horizonDir[0] + Fu[1] * horizonDir[1]);
    var Fv = [Fu[0] + k * horizonDir[0], Fu[1] + k * horizonDir[1]];

    return Fv
}

/**
 * Computes the focal length based on two vanishing points and a center of projection.
 * @see 3.2 "Determining the focal length from a single image"
 * @param Fu the first vanishing point in normalized image coordinates.
 * @param Fv the second vanishing point in normalized image coordinates.
 * @param P the center of projection in normalized image coordinates.
 * @return The relative focal length.
 */
function computeFocalLength(Fu, Fv, P)
{

    //compute Puv, the orthogonal projection of P onto FuFv
    /*
    var dirFuFv = normalize([x - y for x, y in zip(Fu, Fv)])
    var FvP = [x - y for x, y in zip(P, Fv)]
    var proj = dot(dirFuFv, FvP)
    var Puv = [proj * x + y for x, y in zip(dirFuFv, Fv)]

    var PPuv = length([x - y for x, y in zip(P, Puv)])

    var FvPuv = length([x - y for x, y in zip(Fv, Puv)])
    var FuPuv = length([x - y for x, y in zip(Fu, Puv)])
    var FuFv = length([x - y for x, y in zip(Fu, Fv)])*/
    //print("FuFv", FuFv, "FvPuv + FuPuv", FvPuv + FuPuv)
   
    var dirFuFv = normalize(vecSubtract(Fu, Fv));
    var FvP = vecSubtract(Fu, Fv);
    var proj = dot(dirFuFv, FvP);
    var Puv = [proj * dirFuFv[0] + Fv[0], proj * dirFuFv[1] + Fv[1]];

    var PPuv = length(vecSubtract(P, Puv));

    var FvPuv = length(vecSubtract(Fv, Puv));
    var FuPuv = length(vecSubtract(Fu, Puv));
    var FuFv = length(vecSubtract(Fu, Fv));
    

    var fSq = FvPuv * FuPuv - PPuv * PPuv
    //print("FuPuv", FuPuv, "FvPuv", FvPuv, "PPuv", PPuv, "OPuv", FvPuv * FuPuv)
    //print("fSq = ", fSq, " = ", FvPuv * FuPuv, " - ", PPuv * PPuv)
    if (fSq < 0)
    {
        return null;
    }
    
    var f = Math.sqrt(fSq);
    //print("dot 1:", dot(normalize(Fu + [f]), normalize(Fv + [f])))

    return f
}

/**
 *
 */
function computeTriangleOrthocenter(verts)
{
    //print("verts", verts)
    //assert(len(verts) == 3)

    A = verts[0];
    B = verts[1];
    C = verts[2];

    //print("A, B, C", A, B, C)

    a = A[0];
    b = A[1];
    c = B[0];
    d = B[1];
    e = C[0];
    f = C[1];

    N = b * c+ d * e + f * a - c * f - b * e - a * d;
    x = ((d-f) * b * b + (f-b) * d * d + (b-d) * f * f + a * b * (c-e) + c * d * (e-a) + e * f * (a-c) ) / N;
    y = ((e-c) * a * a + (a-e) * c * c + (c-a) * e * e + a * b * (f-d) + c * d * (b-f) + e * f * (d-b) ) / N;
    
    return [x, y];
}

/**
 * Computes the camera rotation matrix based on two vanishing points 
 * and a focal length as in section 3.3 "Computing the rotation matrix".
 * @param Fu the first vanishing point in normalized image coordinates.
 * @param Fv the second vanishing point in normalized image coordinates.
 * @param f the relative focal length.
 * @param The matrix Moc
 */
function computeCameraRotationMatrix(Fu, Fv, f, P)
{

    Fu[0] -= P[0];
    Fu[1] -= P[1];

    Fv[0] -= P[0];
    Fv[1] -= P[1];

    var OFu = [Fu[0], Fu[1], f];
    var OFv = [Fv[0], Fv[1], f];

    //print("matrix dot", dot(OFu, OFv))

    var s1 = length(OFu);
    var upRc = normalize(OFu);

    var s2 = length(OFv);
    var vpRc = normalize(OFv);

    var wpRc = [upRc[1] * vpRc[2] - upRc[2] * vpRc[1], 
                upRc[2] * vpRc[0] - upRc[0] * vpRc[2], 
                upRc[0] * vpRc[1] - upRc[1] * vpRc[0]]; 

    var M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    
    M[0][0] = Fu[0] / s1;
    M[0][1] = Fv[0] / s2;
    M[0][2] = wpRc[0];

    M[1][0] = Fu[1] / s1;
    M[1][1] = Fv[1] / s2;
    M[1][2] = wpRc[1];

    M[2][0] = f / s1;
    M[2][1] = f / s2;
    M[2][2] = wpRc[2];


    return M
}

/**
 * 
 */
function computeVanishingPoint(segment1, segment2)
{
    var x1 = segment1[0][0];
    var y1 = segment1[0][1];
    
    var x2 = segment1[1][0];
    var y2 = segment1[1][1];
    
    var x3 = segment2[0][0];
    var y3 = segment2[0][1];
  
    var x4 = segment2[1][0];
    var y4 = segment2[1][1];
    
    var denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    
    var x = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4);
    x = x / denom;
    
    var y = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4);
    y = y / denom;
    
    return [x, y];
}

function relImgCoords2ImgPlaneCoords(pt, imageWidth, imageHeight)
{
    var ratio = imageWidth / imageHeight; //TODO: float division?
    var sw = ratio;
    var sh = 1;
    return [sw * (pt[0] - 0.5), sh * (pt[1] - 0.5)];
}

function computeCameraParameters()
{
        
        /*def gatherGreasePencilSegments(self):
        '''Collects and returns line segments in normalized image coordinates
        from the first two grease pencil layers.
        \return A list of line segment sets. [i][j][k][l] is coordinate l of point k 
        in segment j from layer i. 
        '''*/
    
        var imageWidth = window.currentImage.naturalWidth;
        var imageHeight = window.currentImage.naturalHeight;
    
        var vpsImPlCoords = [null, null, null];
        for (var i = 0; i < 3; i++)
        {
            var s1 = [controlPoints[4 * i], controlPoints[4 * i + 1]];
            var s2 = [controlPoints[4 * i + 2], controlPoints[4 * i + 3]];
            vanishingPoints[i] = computeVanishingPoint(s1, s2);
            vpsImPlCoords[i] = relImgCoords2ImgPlaneCoords(vanishingPoints[i], imageWidth, imageHeight);
        }
    
        window.computedPrincipalPoint =  computeTriangleOrthocenter(vpsImPlCoords);
        //console.log("window.computedPrincipalPoint " + window.computedPrincipalPoint);

        
        //principal point in image plane coordinates
        var P = window.numVanishingPoints == 3 ? window.computedPrincipalPoint : controlPoints[CP_PRINCIPAL];
        
        if (window.numVanishingPoints == 1)
        {
            //calibration using a single vanishing point
            var imgAspect = imageWidth / imageHeight; //TODO: float div?
            
            //compute the horizon direction
            var horizDir = [1.0, 0.0];
            if (window.useHorizon)
            {
                var xHorizDir = imgAspect * (controlPoints[CP_HORIZ_END][0] - controlPoints[CP_HORIZ_START][0]);
                var yHorizDir = controlPoints[CP_HORIZ_END][1] - controlPoints[CP_HORIZ_START][1];
                horizDir = normalize([-xHorizDir, -yHorizDir])
            }
            //print("horizDir", horizDir)
            
            //compute the vanishing point location
            var vp1 = vanishingPoints[0];
            
            //get the current relative focal length
            /*vfAbs = activeSpace.clip.tracking.camera.focal_length
            sensorWidth = activeSpace.clip.tracking.camera.sensor_width
            
            f = fAbs / sensorWidth * imgAspect
            //print("fAbs", fAbs, "f rel", f)
             */
            var f = 0.7;
            Fu = relImgCoords2ImgPlaneCoords(vp1, imageWidth, imageHeight);
            Fv = computeSecondVanishingPoint(Fu, f, P, horizDir);
        }
        else
        {
            
            
            
            //calibration using two vanishing points
            /*
            if scn.optical_center_type == 'camdata':
                #get the principal point location from camera data
                P = [x for x in  activeSpace.clip.tracking.camera.principal]
                #print("camera data optical center", P[:])
                P[0] /= imageWidth
                P[1] /= imageHeight
                #print("normlz. optical center", P[:])
                P = self.relImgCoords2ImgPlaneCoords(P, imageWidth, imageHeight)
            elif scn.optical_center_type == 'compute':
                if len(vpLineSets) < 3:
                    self.report({'ERROR'}, "A third grease pencil layer is needed to compute the optical center.")
                    return{'CANCELLED'}
                #compute the principal point using a vanishing point from a third gp layer.
                #this computation does not rely on the order of the line sets
                vps = [self.computeIntersectionPointForLineSegments(vpLineSets[i]) for i in range(len(vpLineSets))]
                vps = [self.relImgCoords2ImgPlaneCoords(vps[i], imageWidth, imageHeight) for i in range(len(vps))]
                P = self.computeTriangleOrthocenter(vps)
            else:
                #assume optical center in image midpoint
                pass
            
            #compute the two vanishing points
            vps = [self.computeIntersectionPointForLineSegments(vpLineSets[i]) for i in range(2)]    
        
            #order vanishing points along the image x axis
            if vps[1][0] < vps[0][0]:
                vps.reverse()
                vpLineSets.reverse()
                vpAxisIndices.reverse()            
             */
            //compute focal length
            /*
            
            */
            Fu = relImgCoords2ImgPlaneCoords(vanishingPoints[0], imageWidth, imageHeight);
            Fv = relImgCoords2ImgPlaneCoords(vanishingPoints[1], imageWidth, imageHeight);
            
            var f = computeFocalLength(Fu, Fv, P);
            
        }
        
        //compute camera orientation
        console.log("Fu " + Fu + ", Fv " + Fv + ", f " + f);
        //initial orientation based on the vanishing points and focal length
        //M = computeCameraRotationMatrix(Fu, Fv, f, P)
        
        //sanity check: M should be a pure rotation matrix, 
        //so its determinant should be 1
        /*
        eps = 0.00001
        if 1.0 - M.determinant() < -eps or 1.0 - M.determinant() > eps:
            self.report({'ERROR'}, "Non unit rotation matrix determinant: " + str(M.determinant()))    
            #return{'CANCELLED'} **/
        
        /*
        #align the camera to the coordinate axes as specified
        M = self.alignCoordinateAxes(M, vpAxisIndices[0], vpAxisIndices[1])
        #apply the transform to the camera
        cam.matrix_world = M*/
        
        /*
        move the camera an arbitrary distance away from the ground plane
        TODO: focus on the origin or something
        
        cam.location = (0, 0, 2)*/
    
        /*
        #compute an absolute focal length in mm based 
        #on the current camera settings
        #TODO: make sure this works for all combinations of
        #image dimensions and camera sensor settings
        if imageWidth >= imageHeight:
            fMm = cam.data.sensor_height * f
        else:
            fMm = cam.data.sensor_width * f
        cam.data.lens = fMm
        self.report({'INFO'}, "Camera focal length set to " + str(fMm))
        
        #move principal point of the blender camera
        r = imageWidth / float(imageHeight)
        cam.data.shift_x = -1 * P[0] / r
        cam.data.shift_y = -1 * P[1] / r*/
        
        
        
}

function checkBrowserCompatibility()
{
    var errorList = "";
    
    // Check for the various File API support.
    if (window.File && window.FileReader && window.FileList && window.Blob) 
    {
        //ok
    }
    else 
    {
        errorList += '<li>The File APIs are not fully supported in this browser.</li>';
    }
    
    // check if canvas is supported
    if (window.imageCanvas && window.imageCanvas.getContext &&
        window.overlayCanvas && window.overlayCanvas.getContext) 
    {
        //ok
    }
    else
    {
        errorList += '<li>Canvas not supported.</li>';
    }
    
    return null;
}

function setLabelText(id, text)
{
    $('#' + id).html(text);
}

/**
 *
 */
function startLoadingImageFromURL(url)
{
    var newImage = new Image();
    newImage.src = url;
    
    newImage.onerror = function()
    {
        alert('Error loading image ' + url);
        $('#imageURL').val("");
        newImage = null;
    };
    
    newImage.onload = function()
    {
        window.currentImage = null;
        window.currentImage = newImage;
        draw();
    };
}

/**
 *
 */
function handleFileSelect(evt) 
{
    evt.stopPropagation();
    evt.preventDefault();

    var files = evt.dataTransfer.files; // FileList object.

    if (files.length != 1)
    {
        alert('Try dropping a single image.');
        return;
    }

    // files is a FileList of File objects. List some properties.
    var file = files[0];

    if (typeof FileReader !== "undefined" && file.type.indexOf("image") != -1) 
    {
        window.currentImage = null;
        window.currentImage = new Image();
        var reader = new FileReader();
        reader.onload = function (evt) 
        {
            console.log(evt.target.result);
            window.currentImage.src = evt.target.result;
        };
        reader.readAsDataURL(file);
        
        $('#imageURL').val("")
    }
    else
    {
        alert('Not a recognized image file.');
        return;
    }

    var sc = $("#splashContainer")

    if (sc)
    {
        sc.fadeTo('slow', 0);
    }

    window.currentImage.onload = function()
    {
        draw();
    };
    
    window.currentImage.src = "localfile.png";

}

function controlPointHitTest(x, y)
{
    //hit test in reverse draw order
    for (var i = CP_COUNT - 1; i >= 0; i--)
    {
        var pSc = relImageToScreen(controlPoints[i][0], controlPoints[i][1]);
        var xi = pSc[0];
        var yi = pSc[1];
        
        var distSq = (x - xi) * (x - xi) + (y - yi) * (y - yi);
        
        if (distSq <= SELECTION_RADIUS * SELECTION_RADIUS)
        {
            if (i / 4 >= window.numVanishingPoints &&
                i / 4 < 3)
            {
                //console.log("skipping " + i);
                continue;
            }

            if ((i == CP_HORIZ_END || i == CP_HORIZ_START) &&
                !window.useHorizon)
            {

                continue;
            }

            if (i == CP_PRINCIPAL &&
                window.numVanishingPoints == 3)
            {
                continue;
            }
            
            //console.log("grabbing cp " + i);
            return i;
        }
    }
    
    return -1;
}

/**
 * 
 */
function handleDragOver(evt) 
{
    evt.stopPropagation();
    evt.preventDefault();
    evt.dataTransfer.dropEffect = 'copy'; 
}
    
function onMouseMove(e)
{
    var x = e.pageX - $(this).offset().left;
    var y = e.pageY - $(this).offset().top;
    
    $("#" + ID_OVERLAY_CANVAS).css('cursor', 'normal');
    
    if (window.draggedControlPoint >= 0)
    {
        var idx = window.draggedControlPoint;
        
        var prevDragPosRel = screenToRelImage(window.prevDragPos[0], window.prevDragPos[1]);
        var dragPosRel = screenToRelImage(x, y);
        
        var dx = dragPosRel[0] - prevDragPosRel[0];
        var dy = dragPosRel[1] - prevDragPosRel[1];
        
        window.prevDragPos = [x, y];
        
        controlPoints[idx][0] += dx;
        controlPoints[idx][1] += dy;
        
        controlPoints[idx][0] = Math.min(Math.max(0, controlPoints[idx][0]), 1.0);
        controlPoints[idx][1] = Math.min(Math.max(0, controlPoints[idx][1]), 1.0);
        
        computeCameraParameters();
        
        drawOverlay();
    }
    
    var mouseOverIdx = window.draggedControlPoint >= 0 ? window.draggedControlPoint : controlPointHitTest(x, y);
    
    if (mouseOverIdx >= 0)
    {
        $("#" + ID_OVERLAY_CANVAS).css('cursor', 'pointer');
        var pointName = "";
        if (mouseOverIdx < CP_VLY1_START)
        {
            pointName = "X axis vanishing line";
        }
        else if (mouseOverIdx < CP_VLZ1_START)
        {
            pointName = "Y axis vanishing line";
        }
        else if (mouseOverIdx < CP_HORIZ_START)
        {
            pointName = "Z axis vanishing line";
        }
        else if (mouseOverIdx < CP_PRINCIPAL)
        {
            pointName = "Horizon line";
        }
        else if (mouseOverIdx == CP_PRINCIPAL)
        {
            pointName = "Principal point";
        }
        else if (mouseOverIdx == CP_ORIGIN)
        {
            pointName = "3D origin";
        }
        
        var xHover = truncateFloat(controlPoints[mouseOverIdx][0]);
        var yHover = truncateFloat(controlPoints[mouseOverIdx][1]);
        setLabelText(ID_INFO_LABEL, pointName + " (" + xHover + ", " + yHover + ")");
    }
    else
    {
        setLabelText(ID_INFO_LABEL, "");
    }
    
    if (0)
    {
        var rel = screenToRelImage(x, y);
        var sref = relImageToScreen(rel[0], rel[1]);
        setLabelText(ID_INFO_LABEL, "sc: " + x + ", " + y + " - rel: " + rel[0] + ", " + rel[1] + " - sc ref: " + sref[0] + ", " + sref[1]);
    }
}

function onMouseDown(e)
{
    var x = e.pageX - $(this).offset().left;
    var y = e.pageY - $(this).offset().top; 
    
    window.draggedControlPoint = controlPointHitTest(x, y);
    window.prevDragPos = [x, y];
}

function onMouseEnter(e)
{
    //console.log("Mouse enter");
}

function onMouseLeave(e)
{
    //console.log("Mouse leave");
    window.draggedControlPoint = -1;
}

function onMouseUp(e)
{
    setLabelText(ID_SHARE_URL_LABEL, getShareURL());
    //console.log("Mouse up");
    window.draggedControlPoint = -1;
}

function getCurrentImageRectSc()
{
    var w = $("#" + ID_IMAGE_CANVAS).width();
    var h = $("#" + ID_IMAGE_CANVAS).height();

    var wIm = window.currentImage != null ? window.currentImage.naturalWidth : w;
    var hIm = window.currentImage != null ? window.currentImage.naturalHeight : h;
    
    var wScale = w / wIm;
    var hScale = h / hIm;
    
    var scale = hScale > wScale ? wScale : hScale;
    var wR = scale * wIm;
    var hR = scale * hIm;
    var dx = 0.5 * (w - wR);
    var dy = 0.5 * (h - hR);
    
    return [dx, dy, wR, hR];
}

/**
 * Redraws the background image.
 */
function drawImage()
{
    var w = $("#" + ID_IMAGE_CANVAS).width();
    var h = $("#" + ID_IMAGE_CANVAS).height();
    
    $("#" + ID_IMAGE_CANVAS).attr('width', w);
    $("#" + ID_IMAGE_CANVAS).attr('height', h);
    
    if (window.currentImage)
    {
        var rect = getCurrentImageRectSc();
        var ctx = window.imageCanvas.getContext('2d');
        ctx.drawImage(window.currentImage, rect[0], rect[1], rect[2], rect[3]);
    }
}

/**
 * Redraws the overlay (vanishing lines etc).
 */
function drawOverlay()
{   
    var w = $("#" + ID_OVERLAY_CANVAS).width();
    var h = $("#" + ID_OVERLAY_CANVAS).height();
    
    $("#" + ID_OVERLAY_CANVAS).attr('width', w);
    $("#" + ID_OVERLAY_CANVAS).attr('height', h);
    
    //clear the canvas
    ctx = window.overlayCanvas.getContext('2d');
    ctx.clearRect(0, 0, ctx.canvas.width, ctx.canvas.height);
    ctx.lineWidth = 2;
    
    //draw vanishing lines and vanishing points
    var axisColors = [COL_VLX, COL_VLY, COL_VLZ];
    
    for (var i = 0; i < numVanishingPoints; i++)
    {
        ctx.strokeStyle = axisColors[i];
        ctx.fillStyle = axisColors[i];
        drawLineSegment(ctx, controlPoints[4 * i], controlPoints[4 * i + 1], true);
        drawLineSegment(ctx, controlPoints[4 * i + 2], controlPoints[4 * i + 3], true);
        
        //ctx.lineWidth = 0.5;
        //drawLineSegment(ctx, controlPoints[4 * i], vanishingPoints[i], false);
        //drawLineSegment(ctx, controlPoints[4 * i + 2], vanishingPoints[i], false);
    }
    
    for (var i = 0; i < 3; i++)
    {
        ctx.lineWidth = 0.5;
        ctx.strokeStyle = axisColors[i];
        ctx.fillStyle = axisColors[i];
        drawLineSegment(ctx, controlPoints[CP_ORIGIN], vanishingPoints[i], false);
    }
    
    ctx.lineWidth = 2;
    
    //draw horizon
    if (window.numVanishingPoints == 1 && window.useHorizon)
    {
        ctx.strokeStyle = COL_HORIZ;
        ctx.fillStyle = COL_HORIZ;
        drawLineSegment(ctx, controlPoints[CP_HORIZ_START], controlPoints[CP_HORIZ_END], true);
    }
    
    //draw principal point
    ctx.strokeStyle = COL_PRINCIPAL;
    ctx.fillStyle = COL_PRINCIPAL;
    if (window.numVanishingPoints < 3)
    {
        drawControlPoint(ctx,
                         controlPoints[CP_PRINCIPAL][0],
                         controlPoints[CP_PRINCIPAL][1],
                         true);
    }
    else
    {
        drawControlPoint(ctx,
                         window.computedPrincipalPoint[0],
                         window.computedPrincipalPoint[1],
                         false);
    }
    
    //draw 3D origin
    ctx.strokeStyle = COL_ORIGIN;
    ctx.fillStyle = COL_ORIGIN;
    drawControlPoint(ctx, 
                     controlPoints[CP_ORIGIN][0], 
                     controlPoints[CP_ORIGIN][1], 
                     true);
}

/**
 * 
 */
function drawLineSegment(ctx, start, end, endMarkers)
{
    var p0 = relImageToScreen(start[0], start[1]);
    var p1 = relImageToScreen(end[0], end[1]);
    
    

    ctx.beginPath();
    ctx.moveTo(p0[0], p0[1]);
    ctx.lineTo(p1[0], p1[1]);
    ctx.stroke();
    
    if (endMarkers)
    {
        drawControlPoint(ctx, start[0], start[1], true);
        drawControlPoint(ctx, end[0], end[1], true);
    }
}

/**
 * Draws a single line segment control point at a
 * given position, using a given context.
 */
function drawControlPoint(ctx, x, y, fill)
{
    var p0 = relImageToScreen(x, y);
    ctx.beginPath();
    ctx.arc(p0[0], p0[1], SELECTION_RADIUS, 0 , 2 * Math.PI, false);

    if (fill)
    {
        ctx.fill();
    }
    else
    {
        ctx.stroke();
    }
    
    ctx.closePath();     
}


/**
 * Redraws the image as well as the overlay.
 */
function draw() 
{
    drawImage();
    drawOverlay();
}

function createCameraPresets()
{
    //[Name, sensor width in mm, sensor height in mm]
    window.cameraPresets =
    [
        ["Custom", 36, 24],
        ["APS-C DSLR", 22.3, 14.9],
        ["Canon 1D", 27.9, 18.6],
        ["Canon 1DS", 36, 24],
        ["Canon 5D", 36, 24],
        ["Canon 5D", 36, 24],
        ["Canon 7D", 22.3, 14.9],
        ["Canon 60D", 22.3, 14.9],
        ["Canon 500D", 22.3, 14.9],
        ["Canon 550D", 22.3, 14.9],
        ["Canon 600D", 22.3, 14.9],
        ["Canon 1100D", 22.2, 14.7],
        ["35 mm film", 36, 24],
        ["Micro four thirds", 17.3, 14],
        ["Nikon D3S", 36, 23.9],
        ["Nikon D90", 23.6, 15.8],
        ["Nikon D300S", 23.6, 15.8],
        ["Nikon D3100", 23.1, 15.4],
        ["Nikon D5000", 23.6, 15.8],
        ["Nikon D5100", 23.6, 15.6],
        ["Nikon D7000", 23.6, 15.6],
        ["Red Epic", 30, 15],
        ["Red Epic 2k", 11.1, 6.24],
        ["Red Epic 3k", 16.65, 9.36],
        ["Red Epic 4k", 22.2, 12.6],
        ["Sony A55", 23.4, 15.6],
        ["Super 16 film", 12.52, 7.41],
        ["Super 35 film", 24.89, 18.66]
    ];
    
    //populate the dropdown menu
    var html = '<select id="' + ID_CAMERA_MENU + '">';
    for (var i = 0; i < window.cameraPresets.length; i++)
    {
        html += '<option value="' + i + '">' + window.cameraPresets[i][0] + '</option>';
    }
    html += "</select>";
    
    $("#cameraMenuContainer").html(html);
    
    $('#' + ID_CAMERA_MENU).change(function() 
    {
        setCameraPreset($(this).val());
    });
    
    setCameraPreset(0);
}

function getSensorWidth()
{
    return window.cameraPresets[getCameraPreset()][1];
}

function setCustomSensorWidth(val)
{
    window.cameraPresets[0][1] = val;
}

function getSensorHeight()
{
    return window.cameraPresets[getCameraPreset()][2];
}

function setCustomSensorHeight(val)
{
    window.cameraPresets[0][2] = val;
}

function getCameraPreset()
{
    //console.log("cam preset" + $('#' + ID_CAMERA_MENU).val());
    return $('#' + ID_CAMERA_MENU).val();
}

/**
 *
 */
function setCameraPreset(idx)
{
    console.log("setCameraPreset(idx=" + idx + ")");
    
    $('#' + ID_SENSOR_WIDTH_TEXTFIELD).val(window.cameraPresets[idx][1]);
    $('#' + ID_SENSOR_HEIGHT_TEXTFIELD).val(window.cameraPresets[idx][2]);
    
    if (idx != 0)
    {
        $('#' + ID_SENSOR_WIDTH_TEXTFIELD).attr("disabled", "disabled");
        $('#' + ID_SENSOR_HEIGHT_TEXTFIELD).attr("disabled", "disabled");
    }
    else
    {
        $('#' + ID_SENSOR_WIDTH_TEXTFIELD).removeAttr("disabled");
        $('#' + ID_SENSOR_HEIGHT_TEXTFIELD).removeAttr("disabled");
    }   
}

/**
 *
 */
function setNumVanishingPoints(num)
{
    window.numVanishingPoints = num;
    
    $("input[name='" + ID_NUM_VANISHING_POINTS + '"]')[num - 1].checked = true;
    
    if (num != 1)
    {
        $('#' + ID_HORIZON_CHECKBOX).attr("disabled", "disabled");
    }
    else
    {
        $('#' + ID_HORIZON_CHECKBOX).removeAttr("disabled");
    }
    
    drawOverlay();
}

function getNumVanishingPoints()
{
    return $("input[name='" + ID_NUM_VANISHING_POINTS + "']:checked").val();
}

/**
 *
 */
function setUseHorizonLine(enabled)
{
    window.useHorizon = enabled;
    drawOverlay();
}

/**
 * Returns null on failure, a validated dict on success.
 */
function loadStateFromQueryString(qs)
{
    //construct the dictionary
    var dict = {};
    var query = qs;//window.location.search.substring(1);
    var vars = query.split('&');
    for (var i = 0; i < vars.length; i++)
    {
        var pair = vars[i].split('=');
        var k = pair[0];
        var v = pair[1];
        //console.log(pair[0] + ", " + pair[1]);
        dict[k] = v;
    }
    
    //validate
    var FLOAT = 0;
    var INT = 0;
    var STRING = 1;
    
    var requiredParams =
    [
     PARAM_NUM_VP,
     PARAM_HORIZ_START_X,
     PARAM_HORIZ_START_Y,
     PARAM_HORIZ_END_X,
     PARAM_HORIZ_END_Y,
     PARAM_HORIZ_ENABLED,
     PARAM_CAMERA_PRESET,
     PARAM_SENSOR_WIDTH,
     PARAM_SENSOR_HEIGHT,
     
     PARAM_IMAGE_URL,
     
     PARAM_VLX1_START_X,
     PARAM_VLX1_START_Y,
     PARAM_VLX1_END_X,
     PARAM_VLX1_END_Y,
     
     PARAM_VLX2_START_X,
     PARAM_VLX2_START_Y,
     PARAM_VLX2_END_X,
     PARAM_VLX2_END_Y,
     
     PARAM_VLY1_START_X,
     PARAM_VLY1_START_Y,
     PARAM_VLY1_END_X,
     PARAM_VLY1_END_Y,
     
     PARAM_VLY2_START_X,
     PARAM_VLY2_START_Y,
     PARAM_VLY2_END_X,
     PARAM_VLY2_END_Y,
     
     PARAM_VLZ1_START_X,
     PARAM_VLZ1_START_Y,
     PARAM_VLZ1_END_X,
     PARAM_VLZ1_END_Y,
     
     PARAM_VLZ2_START_X,
     PARAM_VLZ2_START_Y,
     PARAM_VLZ2_END_X,
     PARAM_VLZ2_END_Y,
     
     PARAM_ORIGIN_X,
     PARAM_ORIGIN_Y,
     
     PARAM_PRINCIPAL_X,
     PARAM_PRINCIPAL_Y
     ];
    
    for (var p in requiredParams)
    {
        //TODO: check
    }
    
    //dict is ok. set state
    appState = dict;
    
    //method
    window.numVanishingPoints = appState[PARAM_NUM_VP]
    
    
    setUseHorizonLine(parseInt(appState[PARAM_HORIZ_ENABLED]) == 1);
    
    //camera
    setCameraPreset(appState[PARAM_CAMERA_PRESET]);
    setCustomSensorWidth(appState[PARAM_SENSOR_WIDTH]);
    setCustomSensorHeight(appState[PARAM_SENSOR_HEIGHT]);
    
    if (1)
    {
        controlPoints[CP_HORIZ_START][0] = parseFloat(appState[PARAM_HORIZ_START_X]);
        controlPoints[CP_HORIZ_START][1] = parseFloat(appState[PARAM_HORIZ_START_Y]);
        controlPoints[CP_HORIZ_END][0] = parseFloat(appState[PARAM_HORIZ_END_X]);
        controlPoints[CP_HORIZ_END][1] = parseFloat(appState[PARAM_HORIZ_END_Y]);
        
        //vanishing point positions
        controlPoints[CP_VLX1_START][0] = parseFloat(appState[PARAM_VLX1_START_X]);
        controlPoints[CP_VLX1_START][1] = parseFloat(appState[PARAM_VLX1_START_Y]);
        controlPoints[CP_VLX1_END][0] = parseFloat(appState[PARAM_VLX1_END_X]);
        controlPoints[CP_VLX1_END][1] = parseFloat(appState[PARAM_VLX1_END_Y]);
        
        controlPoints[CP_VLX2_START][0] = parseFloat(appState[PARAM_VLX2_START_X]);
        controlPoints[CP_VLX2_START][1] = parseFloat(appState[PARAM_VLX2_START_Y]);
        controlPoints[CP_VLX2_END][0] = parseFloat(appState[PARAM_VLX2_END_X]);
        controlPoints[CP_VLX2_END][1] = parseFloat(appState[PARAM_VLX2_END_Y]);
        
        controlPoints[CP_VLY1_START][0] = parseFloat(appState[PARAM_VLY1_START_X]);
        controlPoints[CP_VLY1_START][1] = parseFloat(appState[PARAM_VLY1_START_Y]);
        controlPoints[CP_VLY1_END][0] = parseFloat(appState[PARAM_VLY1_END_X]);
        controlPoints[CP_VLY1_END][1] = parseFloat(appState[PARAM_VLY1_END_Y]);
        
        controlPoints[CP_VLY2_START][0] = parseFloat(appState[PARAM_VLY2_START_X]);
        controlPoints[CP_VLY2_START][1] = parseFloat(appState[PARAM_VLY2_START_Y]);
        controlPoints[CP_VLY2_END][0] = parseFloat(appState[PARAM_VLY2_END_X]);
        controlPoints[CP_VLY2_END][1] = parseFloat(appState[PARAM_VLY2_END_Y]);
        
        controlPoints[CP_VLZ1_START][0] = parseFloat(appState[PARAM_VLZ1_START_X]);
        controlPoints[CP_VLZ1_START][1] = parseFloat(appState[PARAM_VLZ1_START_Y]);
        controlPoints[CP_VLZ1_END][0] = parseFloat(appState[PARAM_VLZ1_END_X]);
        controlPoints[CP_VLZ1_END][1] = parseFloat(appState[PARAM_VLZ1_END_Y]);
        
        controlPoints[CP_VLZ2_START][0] = parseFloat(appState[PARAM_VLZ2_START_X]);
        controlPoints[CP_VLZ2_START][1] = parseFloat(appState[PARAM_VLZ2_START_Y]);
        controlPoints[CP_VLZ2_END][0] = parseFloat(appState[PARAM_VLZ2_END_X]);
        controlPoints[CP_VLZ2_END][1] = parseFloat(appState[PARAM_VLZ2_END_Y]);
        
        //origin
        controlPoints[CP_ORIGIN][0] = parseFloat(appState[PARAM_ORIGIN_X]);
        controlPoints[CP_ORIGIN][1] = parseFloat(appState[PARAM_ORIGIN_Y]);
        
        //principal point
        controlPoints[CP_PRINCIPAL][0] = parseFloat(appState[PARAM_PRINCIPAL_X]);
        controlPoints[CP_PRINCIPAL][1] = parseFloat(appState[PARAM_PRINCIPAL_Y]);
    }
    
    //image url
    $('#imageURL').val(appState[PARAM_IMAGE_URL]);
    
    
    startLoadingImageFromURL(appState[PARAM_IMAGE_URL]);
    
}

/**
 * Saves the current state of the app to the 
 * state dictionary.
 */
function saveStateToQueryString()
{
    //method
    appState[PARAM_NUM_VP] = window.numVanishingPoints;
    
    appState[PARAM_HORIZ_START_X] = controlPoints[CP_HORIZ_START][0];
    appState[PARAM_HORIZ_START_Y] = controlPoints[CP_HORIZ_START][1];
    appState[PARAM_HORIZ_END_X] = controlPoints[CP_HORIZ_END][0];
    appState[PARAM_HORIZ_END_Y] = controlPoints[CP_HORIZ_END][1];
    appState[PARAM_HORIZ_ENABLED] = window.useHorizon ? 0 : 1;
    
    //camera
    appState[PARAM_CAMERA_PRESET] = getCameraPreset();
    appState[PARAM_SENSOR_WIDTH] = getSensorWidth();
    appState[PARAM_SENSOR_HEIGHT] = getSensorHeight();
    
    //vanishing point positions
    appState[PARAM_VLX1_START_X] = controlPoints[CP_VLX1_START][0];
    appState[PARAM_VLX1_START_Y] = controlPoints[CP_VLX1_START][1];
    appState[PARAM_VLX1_END_X] = controlPoints[CP_VLX1_END][0];
    appState[PARAM_VLX1_END_Y] = controlPoints[CP_VLX1_END][1];
    
    appState[PARAM_VLX2_START_X] = controlPoints[CP_VLX2_START][0];
    appState[PARAM_VLX2_START_Y] = controlPoints[CP_VLX2_START][1];
    appState[PARAM_VLX2_END_X] = controlPoints[CP_VLX2_END][0];
    appState[PARAM_VLX2_END_Y] = controlPoints[CP_VLX2_END][1];
    
    appState[PARAM_VLY1_START_X] = controlPoints[CP_VLY1_START][0];
    appState[PARAM_VLY1_START_Y] = controlPoints[CP_VLY1_START][1];
    appState[PARAM_VLY1_END_X] = controlPoints[CP_VLY1_END][0];
    appState[PARAM_VLY1_END_Y] = controlPoints[CP_VLY1_END][1];
    
    appState[PARAM_VLY2_START_X] = controlPoints[CP_VLY2_START][0];
    appState[PARAM_VLY2_START_Y] = controlPoints[CP_VLY2_START][1];
    appState[PARAM_VLY2_END_X] = controlPoints[CP_VLY2_END][0];
    appState[PARAM_VLY2_END_Y] = controlPoints[CP_VLY2_END][1];
    
    appState[PARAM_VLZ1_START_X] = controlPoints[CP_VLZ1_START][0];
    appState[PARAM_VLZ1_START_Y] = controlPoints[CP_VLZ1_START][1];
    appState[PARAM_VLZ1_END_X] = controlPoints[CP_VLZ1_END][0];
    appState[PARAM_VLZ1_END_Y] = controlPoints[CP_VLZ1_END][1];
    
    appState[PARAM_VLZ2_START_X] = controlPoints[CP_VLZ2_START][0];
    appState[PARAM_VLZ2_START_Y] = controlPoints[CP_VLZ2_START][1];
    appState[PARAM_VLZ2_END_X] = controlPoints[CP_VLZ2_END][0];
    appState[PARAM_VLZ2_END_Y] = controlPoints[CP_VLZ2_END][1];
    
    //origin
    appState[PARAM_ORIGIN_X] = controlPoints[CP_ORIGIN][0];
    appState[PARAM_ORIGIN_Y] = controlPoints[CP_ORIGIN][1];
    
    //principal
    appState[PARAM_PRINCIPAL_X] = controlPoints[CP_PRINCIPAL][0];
    appState[PARAM_PRINCIPAL_Y] = controlPoints[CP_PRINCIPAL][1];
    
    //image url
    appState[PARAM_IMAGE_URL] = $('#imageURL').val();
    
    //create quert string
    var str = "";
    //console.log(appState);
    for (var k in appState)
    {
        //console.log(k + "=" + appState[k]);
        str = str + k + "=" + appState[k] + "&";
    }
    
    return str.substring(0, str.length - 1);
}

function getShareURL()
{
    return window.location + "?" + saveStateToQueryString();
}

function screenToRelImage(x, y)
{
    var rect = getCurrentImageRectSc();
    
    //position in image coordinates
    var xIm = x - rect[0];
    var yIm = y - rect[1];
    
    var xImRel = Math.min(1.0, Math.max(0.0, xIm / rect[2]));
    var yImRel = Math.min(1.0, Math.max(0.0, (rect[3] - yIm) / rect[3]));
    
    return [xImRel, yImRel];
}

function relImageToScreen(x, y)
{
    var rect = getCurrentImageRectSc();
    //console.log("x, y " + x + ", " + y);
    var xIm = rect[0] + x * rect[2];
    var yIm = rect[1] + (1 - y) * rect[3];
    
    return [xIm, yIm];
}

function reset()
{
    loadStateFromQueryString(DEFAULT_STATE_QUERY_STRING);
    draw();
}

/**
 * Performs setup on document ready.
 */
$(document).ready(function()
{
    //enables/disables various debug output
    window.debug = true;
    
    //used to handle control point dragging
    window.draggedControlPoint = -1;
                  
    window.computedPrincipalPoint = [0.5, 0.5];
    
    //store the three main components
    window.imageCanvas = $("#imageCanvas")[0];
    window.overlayCanvas = $("#overlayCanvas")[0];
    window.splashDiv = $("#splashContainer")[0];
     
    //check HTML5 compatibility
    var errorList = checkBrowserCompatibility();
        
    if (errorList != null)
    {
        //TODO: output errors
        alert(errorList);
        return;
    } 
          
    //trigger redraws on window resize
    $(window).resize(function(e)
    {
        draw();
    });
    
    $('input[name="' + ID_NUM_VANISHING_POINTS + '"]').change(function() 
    {
        setNumVanishingPoints(getNumVanishingPoints());
    });
    
    $('input[name="' + ID_HORIZON_CHECKBOX + '"]').change(function() 
    {
        setUseHorizonLine($(this).is(':checked'));
    });
    
    
    //hook up drag and drop listeners  
    window.splashDiv.addEventListener('dragover', handleDragOver, false);
    window.splashDiv.addEventListener('drop', handleFileSelect, false);
        
    window.overlayCanvas.addEventListener('dragover', handleDragOver, false);
    window.overlayCanvas.addEventListener('drop', handleFileSelect, false);
    
    //other mouse listeners
    $("#overlayCanvas").mousemove(onMouseMove);
    $("#overlayCanvas").mousedown(onMouseDown);
    $("#overlayCanvas").mouseenter(onMouseEnter);
    $("#overlayCanvas").mouseleave(onMouseLeave);
    $("#overlayCanvas").mouseup(onMouseUp);
    
    $("#imageURLSubmit").click(function()
    {
        var url = $('#imageURL').val();
        startLoadingImageFromURL(url);
    });
    
    $("#centerPrincipalPoint").click(function()
    {
        controlPoints[CP_PRINCIPAL][0] = 0.5;
        controlPoints[CP_PRINCIPAL][1] = 0.5;
        drawOverlay();
    });
                  
                  
    $("#reset").click(function()
    {
        reset();
    });
    
    createCameraPresets();
    
    loadStateFromQueryString(DEFAULT_STATE_QUERY_STRING);
    $("#splashContainer").remove();
    draw();
    
    if (window.debug)
    {
        //fade out splash screen immediately if we're in debug mode
        /*$("#splashContainer").fadeTo('slow', 0, function(e)
        {
            
        });*/
        
        startLoadingImageFromURL("webcube.png");
        
        setLabelText(ID_SHARE_URL_LABEL, getShareURL());
        setLabelText(ID_CALIBRATION_INFO_LABEL, "testing info label");
        setLabelText(ID_AXIS_ANGLE_LABEL, "testing axis angle label");
        setLabelText(ID_QUAT_LABEL, "testing quat label");
        setLabelText(ID_EULER_LABEL, "testing euler label");
        //setLabelText(ID_MATRIX_LABEL, "testing matrix label");
        setLabelText(ID_TRANSLATION_LABEL, "testing translation label");
        setLabelText(ID_SHARE_URL_LABEL, "testing url label");
        setLabelText(ID_FIELD_OF_VIEW_LABEL, "testing fov label");
        setLabelText(ID_FOCAL_LENGTH_LABEL, "testing focal length label");
    }
   
});

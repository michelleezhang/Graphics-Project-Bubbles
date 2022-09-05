// Vertex shader program----------------------------------
var VSHADER_SOURCE = 
  'uniform mat4 u_ModelMatrix;\n' +
  'uniform mat4 u_ViewMatrix;\n' + //
  'uniform mat4 u_ProjMatrix;\n' + //
  'attribute vec4 a_Position;\n' +
  'attribute vec4 a_Color;\n' +
  'varying vec4 v_Color;\n' +
  'void main() {\n' +
  '  gl_Position = u_ProjMatrix * u_ViewMatrix * u_ModelMatrix * a_Position;\n' + //
  '  gl_PointSize = 10.0;\n' +
  '  v_Color = a_Color;\n' +
  '}\n';

// Fragment shader program----------------------------------
var FSHADER_SOURCE = 
  'precision mediump float;\n' +
  'varying vec4 v_Color;\n' +
  'void main() {\n' +
  '  gl_FragColor = v_Color;\n' +
  '}\n';


// Global Variables
//------------For WebGL-----------------------------------------------
var gl; 
var canvas = document.getElementById('webgl');

// ---------- For shapes & their matrices ---------------------------------
var floatsPerVertex = 7;	// # of Float32Array elements used for each vertex

var modelMatrix = new Matrix4();

var floatsPerVertex = 7;	// # of Float32Array elements used for each vertex
													// (x,y,z,w)position + (r,g,b)color
													// Later, see if you can add:
													// (x,y,z) surface normal + (tx,ty) texture addr.

//------------For Animation---------------------------------------------
var g_isRun = true;                 // run/stop for animation; used in tick().
var g_lastMS = Date.now();    			// Timestamp for most-recently-drawn image; 
                                    // in milliseconds; used by 'animate()' fcn 
                                    // (now called 'timerAll()' ) to find time
                                    // elapsed since last on-screen image.
var g_angle01 = 0;                  // initial rotation angle
var g_angle01Rate = 30.0;           // rotation speed, in degrees/second 

// All of our time-dependent params
var g_angle0now  =   0;       // init Current rotation angle, in degrees
var g_angle0rate = -30.0;       // init Rotation angle rate, in degrees/second.
var g_angle0brake=	 1.0;				// init Speed control; 0=stop, 1=full speed.
var g_angle0min  = -20.0;       // init min, max allowed angle, in degrees.
var g_angle0max  =  20.0;
                                //---------------
var g_angle1now  =   0.0; 		
var g_angle1rate =  30.0;		
var g_angle1brake=	 1.0;			
var g_angle1min  = -20.0;     
var g_angle1max  =  20.0;

var g_angle2now  =   0.0; 
var g_angle2rate =  40.0;
var g_angle2brake=	 1.0;	
var g_angle2min  = -20.0;  
var g_angle2max  = 20.0;

var g_angle3now  =   0.0;
var g_angle3rate =  0.5;	
var g_angle3brake=	 1.0;	
var g_angle3min  = -0.8;       
var g_angle3max  = 0.8;			

//------------For mouse click-and-drag: -------------------------------
var isDrag=false;		// mouse-drag: true when user holds down mouse button
var xMclik=0.0;			// last mouse button-down position (in CVV coords)
var yMclik=0.0;   
var xMdragTot=0.0;	// total (accumulated) mouse-drag amounts (in CVV coords).
var yMdragTot=0.0;  

var qNew = new Quaternion(0,0,0,1); // most-recent mouse drag's rotation
var qTot = new Quaternion(0,0,0,1);	// 'current' orientation (made from qNew)
var quatMatrix = new Matrix4();				// rotation matrix, made from latest qTot



var eyeX = 5.0, eyeY = 0.0, eyeZ = 1.0; // Eye position
var thetarad = Math.PI; //radians
var deltaz = 0.0;
var upX = 0.0, upY = 0.0, upZ = 1.0;


function main() {
//==============================================================================
  // Retrieve <canvas> element
  var canvas = document.getElementById('webgl');

  // Get the rendering context for WebGL
  var gl = getWebGLContext(canvas);
  if (!gl) {
    console.log('Failed to get the rendering context for WebGL');
    return;
  }

  // Initialize shaders
  if (!initShaders(gl, VSHADER_SOURCE, FSHADER_SOURCE)) {
    console.log('Failed to intialize shaders.');
    return;
  }

  // Initialize a Vertex Buffer in the graphics system to hold our vertices
  var n = initVertexBuffer(gl);
  if (n < 0) {
    console.log('Failed to set the vertex information');
    return;
  }

  // Specify the color for clearing <canvas>
  gl.clearColor(0.0, 0.5, 0.5, 0.7);

	// Enable 3D depth-test when drawing: don't over-draw at any pixel 
	// unless the new Z value is closer to the eye than the old one..
	gl.enable(gl.DEPTH_TEST); 	 
	
  // Get handle to graphics system's storage location of u_ModelMatrix

  var u_ProjMatrix = gl.getUniformLocation(gl.program, 'u_ProjMatrix'); // NEW
  var u_ViewMatrix = gl.getUniformLocation(gl.program, 'u_ViewMatrix'); // NEW

  var u_ModelMatrix = gl.getUniformLocation(gl.program, 'u_ModelMatrix');
  if (!u_ModelMatrix || !u_ProjMatrix || !u_ViewMatrix) { 
    console.log('Failed to get the storage location of u_ModelMatrix');
    return;
  }


  // Create a local version of our model matrix in JavaScript 
  var modelMatrix = new Matrix4();
  var viewMatrix = new Matrix4();
  var projMatrix = new Matrix4();

  	// Register the Mouse & Keyboard Event-handlers-------------------------------
	canvas.onmousedown	=	function(ev){myMouseDown( ev, gl, canvas) }; 
	// when user's mouse button goes down, call mouseDown() function
	canvas.onmousemove = 	function(ev){myMouseMove( ev, gl, canvas) };
						  // when the mouse moves, call mouseMove() function					
	canvas.onmouseup = 		function(ev){myMouseUp(   ev, gl, canvas)};
	document.onkeydown = function(ev){ keydown(ev, gl); };

//-----------------  
  // Start drawing: create 'tick' variable whose value is this function:
  var tick = function() {
	requestAnimationFrame(tick, canvas);   
	timerAll();

	resize();
  };
  tick();							// start (and continue) animation: draw current image
	
	function resize() {
		console.log('Canvas width,height=', canvas.width, canvas.height);		
		console.log('Browser window: innerWidth,innerHeight=', innerWidth, innerHeight);	
		
		//Make canvas fill the top 2/3 of our browser window:
		var xtraMargin = 16;    // keep a margin (otherwise, browser adds scroll-bars)
		canvas.width = innerWidth - xtraMargin;
		canvas.height = (innerHeight* (2/3)) - xtraMargin;

		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

		var vpAspect = canvas.width / canvas.height;

		gl.viewport(0, 0, canvas.width/2, canvas.height);

		var near = 1.0;
		var far = 15.0; // 10
		var height = Math.tan((35/2)*(Math.PI)/180) * ((far - near)/3);
		var width = vpAspect * height;

		projMatrix.setPerspective(35.0, // field of view angle
			vpAspect, // aspect ratio
			near, // distance to the near plane
			far); // distance to the far plane

		// calculate view matrix -- eye point, lookat point, up vector
		viewMatrix.setLookAt(eyeX, eyeY, eyeZ, 
			eyeX + Math.cos(thetarad), eyeY + Math.sin(thetarad), eyeZ + deltaz,
			upX, upY, upZ);
		drawAll(gl, n, modelMatrix, u_ModelMatrix, viewMatrix, u_ViewMatrix, projMatrix, u_ProjMatrix);

		gl.viewport(canvas.width/2, 0, canvas.width/2, canvas.height);

		projMatrix.setOrtho(-width, width,
			-height, height,
			near, far); 
		drawAll(gl, n, modelMatrix, u_ModelMatrix, viewMatrix, u_ViewMatrix, projMatrix, u_ProjMatrix);
	}
}


function timerAll() {
	  var nowMS = Date.now();             
	  var elapsedMS = nowMS - g_lastMS;  
	  g_lastMS = nowMS;           
	  if(elapsedMS > 1000.0) {            
		elapsedMS = 1000.0/30.0;
		}
	  g_angle01 = animate(g_angle01);
	  // Find new time-dependent parameters using the current or elapsed time:
	  g_angle0now += g_angle0rate * g_angle0brake * (elapsedMS * 0.001);	// update.
	  g_angle1now += g_angle1rate * g_angle1brake * (elapsedMS * 0.001);
	  g_angle2now += g_angle2rate * g_angle2brake * (elapsedMS * 0.001);
	  g_angle3now += g_angle3rate * g_angle3brake * (elapsedMS * 0.001);
	  // apply angle limits:  going above max, or below min? reverse direction!
	  // (!CAUTION! if max < min, then these limits do nothing...)
	  if((g_angle0now >= g_angle0max && g_angle0rate > 0) || // going over max, or
		   (g_angle0now <= g_angle0min && g_angle0rate < 0)  ) // going under min ?
		   g_angle0rate *= -1;	// YES: reverse direction.
	  if((g_angle1now >= g_angle1max && g_angle1rate > 0) || // going over max, or 
		   (g_angle1now <= g_angle1min && g_angle1rate < 0) )	 // going under min ?
		   g_angle1rate *= -1;	// YES: reverse direction.
	  if((g_angle2now >= g_angle2max && g_angle2rate > 0) || // going over max, or
		   (g_angle2now <= g_angle2min && g_angle2rate < 0) )	 // going under min ?
		   g_angle2rate *= -1;	// YES: reverse direction.
	  if((g_angle3now >= g_angle3max && g_angle3rate > 0) || // going over max, or
		   (g_angle3now <= g_angle3min && g_angle3rate < 0) )	 // going under min ?
		   g_angle3rate *= -1;	// YES: reverse direction.
		// *NO* limits? Don't let angles go to infinity! cycle within -180 to +180.
		if(g_angle0min > g_angle0max)	
		{// if min and max don't limit the angle, then
			if(     g_angle0now < -180.0) g_angle0now += 360.0;	// go to >= -180.0 or
			else if(g_angle0now >  180.0) g_angle0now -= 360.0;	// go to <= +180.0
		}
		if(g_angle1min > g_angle1max)
		{
			if(     g_angle1now < -180.0) g_angle1now += 360.0;	// go to >= -180.0 or
			else if(g_angle1now >  180.0) g_angle1now -= 360.0;	// go to <= +180.0
		}
	  if(g_angle2min > g_angle2max)
		{
			if(     g_angle2now < -180.0) g_angle2now += 360.0;	// go to >= -180.0 or
			else if(g_angle2now >  180.0) g_angle2now -= 360.0;	// go to <= +180.0
		}
	  if(g_angle3min > g_angle3max)
		{
			if(     g_angle3now < -180.0) g_angle3now += 360.0;	// go to >= -180.0 or
			else if(g_angle3now >  180.0) g_angle3now -= 360.0;	// go to <= +180.0
		}
	 }

function initVertexBuffer(gl) {
//==============================================================================
// Create one giant vertex buffer object (VBO) that holds all vertices for all
// shapes.
 
 	// Make each 3D shape in its own array of vertices:
  makeSphere();						// create, fill the sphVerts array
  makeTorus();						// create, fill the torVerts array
  makeTetrahedron();
  makeFishBody();
  makeFishWing();
  makeAxes();
  makeGroundGrid();		

  // how many floats total needed to store all shapes?
	var mySiz = (sphVerts.length + torVerts.length + tetraVerts.length 
		+ fishbodVerts.length + fishwingVerts.length + axesVerts.length + gndVerts.length);
							 	

	// How many vertices total?
	var nn = mySiz / floatsPerVertex;
	console.log('nn is', nn, 'mySiz is', mySiz, 'floatsPerVertex is', floatsPerVertex);
	// Copy all shapes into one big Float32 array:
	
  var colorShapes = new Float32Array(mySiz);
	// Copy them:  remember where to start for each shape:
		sphStart = 0;						// next, we'll store the sphere;
	for(i=0, j=0; j< sphVerts.length; i++, j++) {// don't initialize i -- reuse it!
		colorShapes[i] = sphVerts[j];
		}

		torStart = i;						// next, we'll store the torus;
	for(j=0; j< torVerts.length; i++, j++) {
		colorShapes[i] = torVerts[j];
		}

		tetraStart = i;						
	for(j=0; j< tetraVerts.length; i++, j++) {
		colorShapes[i] = tetraVerts[j];
		}

		fishbodStart = i;						
	for(j=0; j< fishbodVerts.length; i++, j++) {
		colorShapes[i] = fishbodVerts[j];
		}

		fishwingStart = i;						
	for(j=0; j< fishwingVerts.length; i++, j++) {
		colorShapes[i] = fishwingVerts[j];
		}

		axesStart = i;						
	for(j=0; j< axesVerts.length; i++, j++) {
		colorShapes[i] = axesVerts[j];
		}

		gndStart = i;						// next we'll store the ground-plane;
	for(j=0; j< gndVerts.length; i++, j++) {
		colorShapes[i] = gndVerts[j];
		}
  // Create a buffer object on the graphics hardware:
  var shapeBufferHandle = gl.createBuffer();  
  if (!shapeBufferHandle) {
    console.log('Failed to create the shape buffer object');
    return false;
  }

  // Bind the the buffer object to target:
  gl.bindBuffer(gl.ARRAY_BUFFER, shapeBufferHandle);
  // Transfer data from Javascript array colorShapes to Graphics system VBO
  gl.bufferData(gl.ARRAY_BUFFER, colorShapes, gl.STATIC_DRAW);
    
  // Get graphics system's handle for our Vertex Shader's position-input variable: 
  var a_Position = gl.getAttribLocation(gl.program, 'a_Position');
  if (a_Position < 0) {
    console.log('Failed to get the storage location of a_Position');
    return -1;
  }

  var FSIZE = colorShapes.BYTES_PER_ELEMENT; // how many bytes per stored value?

  // Use handle to specify how to retrieve **POSITION** data from our VBO:
  gl.vertexAttribPointer(
  		a_Position, 	// choose Vertex Shader attribute to fill with data
  		4, 						// how many values? 1,2,3 or 4.  (we're using x,y,z,w)
  		gl.FLOAT, 		// data type for each value: usually gl.FLOAT
  		false, 				// did we supply fixed-point data AND it needs normalizing?
  		FSIZE * floatsPerVertex, // Stride -- how many bytes used to store each vertex?
  									// (x,y,z,w, r,g,b) * bytes/value
  		0);						// Offset -- now many bytes from START of buffer to the
  									// value we will actually use?
  gl.enableVertexAttribArray(a_Position);  
  									// Enable assignment of vertex buffer object's position data

  // Get graphics system's handle for our Vertex Shader's color-input variable;
  var a_Color = gl.getAttribLocation(gl.program, 'a_Color');
  if(a_Color < 0) {
    console.log('Failed to get the storage location of a_Color');
    return -1;
  }
  // Use handle to specify how to retrieve **COLOR** data from our VBO:
  gl.vertexAttribPointer(
  	a_Color, 				// choose Vertex Shader attribute to fill with data
  	3, 							// how many values? 1,2,3 or 4. (we're using R,G,B)
  	gl.FLOAT, 			// data type for each value: usually gl.FLOAT
  	false, 					// did we supply fixed-point data AND it needs normalizing?
  	FSIZE * 7, 			// Stride -- how many bytes used to store each vertex?
  									// (x,y,z,w, r,g,b) * bytes/value
  	FSIZE * 4);			// Offset -- how many bytes from START of buffer to the
  									// value we will actually use?  Need to skip over x,y,z,w
  									
  gl.enableVertexAttribArray(a_Color);  
  									// Enable assignment of vertex buffer object's position data

  // Unbind the buffer object 
  gl.bindBuffer(gl.ARRAY_BUFFER, null);

  return nn;
}

function makeSphere() {
//==============================================================================
// Make ring-like equal-lattitude 'slices' of the sphere (bounded by planes of constant z), 
// and connect them to build the sphere from one triangle strip.
  var slices = 13;			// # of slices of the sphere along the z axis. >=3 req'd
											// (choose odd # or prime# to avoid accidental symmetry)
  var sliceVerts	= 27;	// # of vertices around the top edge of the slice
											// (same number of vertices on bottom of slice, too)
  var topColr = new Float32Array([0.7, 0.7, 0.7]);	// North Pole: light gray
  var equColr = new Float32Array([0.3, 0.7, 0.3]);	// Equator:    bright green
  var botColr = new Float32Array([0.9, 0.9, 0.9]);	// South Pole: brightest gray.
  var sliceAngle = Math.PI/slices;	// lattitude angle spanned by one slice.

	// Create a (global) array to hold this sphere's vertices:
  sphVerts = new Float32Array(  ((slices * 2* sliceVerts) -2) * floatsPerVertex);
										// # of vertices * # of elements needed to store them. 
										// each slice requires 2*sliceVerts vertices except 1st and
										// last ones, which require only 2*sliceVerts-1.
										
	// Create dome-shaped top slice of sphere at z=+1
	// s counts slices; v counts vertices; 
	// j counts array elements (vertices * elements per vertex)
	var cos0 = 0.0;					// sines,cosines of slice's top, bottom edge.
	var sin0 = 0.0;
	var cos1 = 0.0;
	var sin1 = 0.0;	
	var j = 0;							// initialize our array index
	var isLast = 0;
	var isFirst = 1;
	for(s=0; s<slices; s++) {	// for each slice of the sphere,
		// find sines & cosines for top and bottom of this slice
		if(s==0) {
			isFirst = 1;	// skip 1st vertex of 1st slice.
			cos0 = 1.0; 	// initialize: start at north pole.
			sin0 = 0.0;
		}
		else {					// otherwise, new top edge == old bottom edge
			isFirst = 0;	
			cos0 = cos1;
			sin0 = sin1;
		}								// & compute sine,cosine for new bottom edge.
		cos1 = Math.cos((s+1)*sliceAngle);
		sin1 = Math.sin((s+1)*sliceAngle);
		// go around the entire slice, generating TRIANGLE_STRIP verts
		// (Note we don't initialize j; grows with each new attrib,vertex, and slice)
		if(s==slices-1) isLast=1;	// skip last vertex of last slice.
		for(v=isFirst; v< 2*sliceVerts-isLast; v++, j+=floatsPerVertex) {	
			if(v%2==0)
			{				// put even# vertices at the the slice's top edge
							// (why PI and not 2*PI? because 0 <= v < 2*sliceVerts
							// and thus we can simplify cos(2*PI(v/2*sliceVerts))  
				sphVerts[j  ] = sin0 * Math.cos(Math.PI*(v)/sliceVerts); 	
				sphVerts[j+1] = sin0 * Math.sin(Math.PI*(v)/sliceVerts);	
				sphVerts[j+2] = cos0;		
				sphVerts[j+3] = 1.0;			
			}
			else { 	// put odd# vertices around the slice's lower edge;
							// x,y,z,w == cos(theta),sin(theta), 1.0, 1.0
							// 					theta = 2*PI*((v-1)/2)/capVerts = PI*(v-1)/capVerts
				sphVerts[j  ] = sin1 * Math.cos(Math.PI*(v-1)/sliceVerts);		// x
				sphVerts[j+1] = sin1 * Math.sin(Math.PI*(v-1)/sliceVerts);		// y
				sphVerts[j+2] = cos1;																				// z
				sphVerts[j+3] = 1.0;																				// w.		
			}
			if(s==0) {	// finally, set some interesting colors for vertices:
				sphVerts[j+4]=topColr[0]; 
				sphVerts[j+5]=topColr[1]; 
				sphVerts[j+6]=topColr[2];	
				}
			else if(s==slices-1) {
				sphVerts[j+4]=botColr[0]; 
				sphVerts[j+5]=botColr[1]; 
				sphVerts[j+6]=botColr[2];	
			}
			else {
					sphVerts[j+4]=Math.random();// equColr[0]; 
					sphVerts[j+5]=Math.random();// equColr[1]; 
					sphVerts[j+6]=Math.random();// equColr[2];					
			}
		}
	}
}

function makeTorus() {
var rbend = 1.0;										// Radius of circle formed by torus' bent bar
var rbar = 0.5;											// radius of the bar we bent to form torus
var barSlices = 23;									// # of bar-segments in the torus: >=3 req'd;
																		// more segments for more-circular torus
var barSides = 13;										// # of sides of the bar (and thus the 
																		// number of vertices in its cross-section)
																		// >=3 req'd;
																		// more sides for more-circular cross-section

	// Create a (global) array to hold this torus's vertices:
 torVerts = new Float32Array(floatsPerVertex*(2*barSides*barSlices +2));
//	Each slice requires 2*barSides vertices, but 1st slice will skip its first 
// triangle and last slice will skip its last triangle. To 'close' the torus,
// repeat the first 2 vertices at the end of the triangle-strip.  Assume 7

var phi=0, theta=0;										// begin torus at angles 0,0
var thetaStep = 2*Math.PI/barSlices;	// theta angle between each bar segment
var phiHalfStep = Math.PI/barSides;		// half-phi angle between each side of bar
																			// (WHY HALF? 2 vertices per step in phi)
	// s counts slices of the bar; v counts vertices within one slice; j counts
	// array elements (Float32) (vertices*#attribs/vertex) put in torVerts array.
	for(s=0,j=0; s<barSlices; s++) {		// for each 'slice' or 'ring' of the torus:
		for(v=0; v< 2*barSides; v++, j+=7) {		// for each vertex in this slice:
			if(v%2==0)	{	// even #'d vertices at bottom of slice,
				torVerts[j  ] = (rbend + rbar*Math.cos((v)*phiHalfStep)) * 
																						 Math.cos((s)*thetaStep);
							  //	x = (rbend + rbar*cos(phi)) * cos(theta)
				torVerts[j+1] = (rbend + rbar*Math.cos((v)*phiHalfStep)) *
																						 Math.sin((s)*thetaStep);
								//  y = (rbend + rbar*cos(phi)) * sin(theta) 
				torVerts[j+2] = -rbar*Math.sin((v)*phiHalfStep);
								//  z = -rbar  *   sin(phi)
				torVerts[j+3] = 1.0;		// w
			}
			else {				// odd #'d vertices at top of slice (s+1);
										// at same phi used at bottom of slice (v-1)
				torVerts[j  ] = (rbend + rbar*Math.cos((v-1)*phiHalfStep)) * 
																						 Math.cos((s+1)*thetaStep);
							  //	x = (rbend + rbar*cos(phi)) * cos(theta)
				torVerts[j+1] = (rbend + rbar*Math.cos((v-1)*phiHalfStep)) *
																						 Math.sin((s+1)*thetaStep);
								//  y = (rbend + rbar*cos(phi)) * sin(theta) 
				torVerts[j+2] = -rbar*Math.sin((v-1)*phiHalfStep);
								//  z = -rbar  *   sin(phi)
				torVerts[j+3] = 1.0;		// w
			}
			torVerts[j+4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j+5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j+6] = Math.random();		// random color 0.0 <= B < 1.0
		}
	}
	// Repeat the 1st 2 vertices of the triangle strip to complete the torus:
			torVerts[j  ] = rbend + rbar;	// copy vertex zero;
						  //	x = (rbend + rbar*cos(phi==0)) * cos(theta==0)
			torVerts[j+1] = 0.0;
							//  y = (rbend + rbar*cos(phi==0)) * sin(theta==0) 
			torVerts[j+2] = 0.0;
							//  z = -rbar  *   sin(phi==0)
			torVerts[j+3] = 1.0;		// w
			torVerts[j+4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j+5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j+6] = Math.random();		// random color 0.0 <= B < 1.0
			j+=7; // go to next vertex:
			torVerts[j  ] = (rbend + rbar) * Math.cos(thetaStep);
						  //	x = (rbend + rbar*cos(phi==0)) * cos(theta==thetaStep)
			torVerts[j+1] = (rbend + rbar) * Math.sin(thetaStep);
							//  y = (rbend + rbar*cos(phi==0)) * sin(theta==thetaStep) 
			torVerts[j+2] = 0.0;
							//  z = -rbar  *   sin(phi==0)
			torVerts[j+3] = 1.0;		// w
			torVerts[j+4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j+5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j+6] = Math.random();		// random color 0.0 <= B < 1.0
}

function makeGroundGrid() {
//==============================================================================
// Create a list of vertices that create a large grid of lines in the x,y plane
// centered at x=y=z=0.  Draw this shape using the GL_LINES primitive.

	var xcount = 100;			// # of lines to draw in x,y to make the grid.
	var ycount = 100;		
	var xymax	= 100.0;			// grid size; extends to cover +/-xymax in x and y.
 	var xColr = new Float32Array([1.0, 1.0, 0.3]);	// bright yellow
 	var yColr = new Float32Array([0.5, 1.0, 0.5]);	// bright green.
 	
	// Create an (global) array to hold this ground-plane's vertices:
	gndVerts = new Float32Array(floatsPerVertex*2*(xcount+ycount));
						// draw a grid made of xcount+ycount lines; 2 vertices per line.
						
	var xgap = xymax/(xcount-1);		// HALF-spacing between lines in x,y;
	var ygap = xymax/(ycount-1);		// (why half? because v==(0line number/2))
	
	// First, step thru x values as we make vertical lines of constant-x:
	for(v=0, j=0; v<2*xcount; v++, j+= floatsPerVertex) {
		if(v%2==0) {	// put even-numbered vertices at (xnow, -xymax, 0)
			gndVerts[j  ] = -xymax + (v  )*xgap;	// x
			gndVerts[j+1] = -xymax;								// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		else {				// put odd-numbered vertices at (xnow, +xymax, 0).
			gndVerts[j  ] = -xymax + (v-1)*xgap;	// x
			gndVerts[j+1] = xymax;								// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		gndVerts[j+4] = xColr[0];			// red
		gndVerts[j+5] = xColr[1];			// grn
		gndVerts[j+6] = xColr[2];			// blu
	}
	// Second, step thru y values as wqe make horizontal lines of constant-y:
	// (don't re-initialize j--we're adding more vertices to the array)
	for(v=0; v<2*ycount; v++, j+= floatsPerVertex) {
		if(v%2==0) {		// put even-numbered vertices at (-xymax, ynow, 0)
			gndVerts[j  ] = -xymax;								// x
			gndVerts[j+1] = -xymax + (v  )*ygap;	// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		else {					// put odd-numbered vertices at (+xymax, ynow, 0).
			gndVerts[j  ] = xymax;								// x
			gndVerts[j+1] = -xymax + (v-1)*ygap;	// y
			gndVerts[j+2] = 0.0;									// z
			gndVerts[j+3] = 1.0;									// w.
		}
		gndVerts[j+4] = yColr[0];			// red
		gndVerts[j+5] = yColr[1];			// grn
		gndVerts[j+6] = yColr[2];			// blu
	}
}

function makeTetrahedron() {
	var c30 = Math.sqrt(0.75);	
	var sq2	= Math.sqrt(2.0);	
	tetraVerts = new Float32Array([
		 0.0,	 0.0, sq2, 1.0,		  1.0, 0.5, 0.9,	// Node 0 light PINK
		 c30, -0.5, 0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 1 PINK
		 0.0,  1.0, 0.0, 1.0,  		0.3, 	0.9,	0.5,	// Node 2 GREEN
		   // Face 1: (right side)
		   0.0,	 0.0, sq2, 1.0,			1.0, 0.5, 0.9,	// Node 0
		 0.0,  1.0, 0.0, 1.0,  		0.3, 	0.9,	0.5,	// Node 2
		 -c30, -0.5, 0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 3
			// Face 2: (lower side)
			 0.0,	 0.0, sq2, 1.0,			1.0, 0.5, 0.9,	// Node 0 
		 -c30, -0.5, 0.0, 1.0, 		1.0, 0.5, 0.9, 	// Node 3
		 c30, -0.5, 0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 1 
			 // Face 3: (base side)  
		 -c30, -0.5,  0.0, 1.0, 		1.0, 0.5, 0.9, 	// Node 3
		 0.0,  1.0,  0.0, 1.0,  	0.3, 	0.9,	0.5,	// Node 2
		 c30, -0.5,  0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 1
	
		 //BOTTOM
		 -c30, -0.5, 0.0, 1.0, 		1.0, 0.5, 0.9, 	// Node 3
		 c30, -0.5,  0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 1
		 0.0,	 0.0, -sq2, 1.0,			0.3, 	0.9,	0.5,	// Node 4 which is rlly just node 0 shifted down
	
		 0.0,	 0.0, -sq2, 1.0,		0.5, 0.5, 0.9,	// Node 4
		 0.0,  1.0,  0.0, 1.0,  	0.3, 	0.9,	0.5,	// Node 2
		 -c30, -0.5, 0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 3
	
		 c30, -0.5,  0.0, 1.0, 		0.8, 	0.09,	0.6, 	// Node 1 	0.8, 	0.09,	0.6,
		 0.0,  1.0,  0.0, 1.0,  	0.3, 	0.9,	0.5,	// Node 2 0.3, 	0.9,	0.5,
		 0.0,	 0.0, -sq2, 1.0,		0.5, 0.5, 0.9,	// Node 4 1.0, 0.5, 0.9,
		]);
}

function makeFishBody() {
	  fishbodVerts = new Float32Array([
	  // HEXAGON PRISM
      // bottom base
	  0.0, 1.0, 0.5, 1.0,			0.2, 0.1, 0.0,  // point NODE 1 WHITE
	  -0.5, 0.5, 0.0, 1.0,			1.0,	0.6,	0.0, // node 2 ORANGE
	  0.5, 0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 3 PURPLE
   
	  0.5, 0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 3
	  -0.5, 0.5, 0.0, 1.0,			1.0,	0.6,	0.0, // node 2
	  0.5, -0.5, 0.0, 1.0,			0.2, 0.1, 0.0, // node 4
   
	  0.5, -0.5, 0.0, 1.0,		0.2, 0.1, 0.0, // node 4 ORANGE
	  -0.5, 0.5, 0.0, 1.0,			1.0,	0.6,	0.0, // node 2
	  -0.5, -0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 5
   
	  -0.5, -0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 5
	  0.5, -0.5, 0.0, 1.0,			0.2, 0.1, 0.0, // node 4
	  0.0, -1.0, 0.5, 1.0,			1.0,	0.6,	0.0, //point NODE 6
	  
	  // top base
	  0.0, 1.0, 0.5, 1.0,			0.2, 0.1, 0.0,  // NODE 1
	  -0.5, 0.5, 1.0, 1.0,			1.0, 	0.6,	0.0, // node 2 MOVED UP
	  0.5, 0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 3
   
	  0.5, 0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 3
	  -0.5, 0.5, 1.0, 1.0,			1.0, 	0.6,	0.0, // node 2
	  0.5, -0.5, 1.0, 1.0,			0.2, 0.1, 0.0, // node 4
   
	  0.5, -0.5, 1.0, 1.0,			0.2, 0.1, 0.0, // node 4
	  -0.5, 0.5, 1.0, 1.0,			1.0, 	0.6,	0.0, // node 2
	  -0.5, -0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 5
   
	  -0.5, -0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 5
	  0.5, -0.5, 1.0, 1.0,			0.2, 0.1, 0.0, // node 4
	  0.0, -1.0, 0.5, 1.0,			1.0, 	0.6,	0.0, // node 6
   
	  // side base 1
	  0.0, 1.0, 0.5, 1.0,			0.2, 0.1, 0.0,  // node 1
	  -0.5, 0.5, 0.0, 1.0,			1.0, 	0.6,	0.0, // node 2 B
	  -0.5, 0.5, 1.0, 1.0,			1.0, 	0.6,	0.0, //DONE node 2 T
   
	  -0.5, 0.5, 1.0, 1.0,			1.0, 	0.6,	0.0, // node 2 T
	  -0.5, 0.5, 0.0, 1.0,			1.0, 	0.6,	0.0, // node 2 B
	  -0.5, -0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 5 T
   
	  -0.5, -0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 5 T
	  -0.5, 0.5, 0.0, 1.0,			1.0, 	0.6,	0.0, // node 2 B
	  -0.5, -0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 5 B
   
	  -0.5, -0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 5 B
	  -0.5, -0.5, 1.0, 1.0,			0.6, 0.7, 0.9, // node 5 T
	  0.0, -1.0, 0.5, 1.0,			1.0, 	0.6,	0.0, // node 6
   
	  // side base 2
	  0.0, 1.0, 0.5, 1.0,			0.2, 0.1, 0.0,  // node 1 
	  0.5, 0.5, 0.0, 1.0,			0.6, 0.7, 0.9, // node 3 B
	  0.5, 0.5, 1.0, 1.0,		 0.6, 0.7, 0.9, // node 3 T
   
	  0.5, 0.5, 1.0, 1.0,        0.6, 0.7, 0.9, // node 3 T
	  0.5, 0.5, 0.0, 1.0,        0.6, 0.7, 0.9, // node 3 B
	  0.5, -0.5, 1.0, 1.0,       1.0, 	0.6,	0.0, // node 4 T
	  
	  0.5, -0.5, 1.0, 1.0,       1.0, 	0.6,	0.0,  // node 4 T
	  0.5, 0.5, 0.0, 1.0,        0.6, 0.7, 0.9, // node 3 B
	  0.5, -0.5, 0.0, 1.0,			1.0, 	0.6,	0.0, // node 4 B
   
	  0.5, -0.5, 0.0, 1.0,			1.0, 	0.6,	0.0, // node 4 B
	  0.5, -0.5, 1.0, 1.0,			1.0, 	0.6,	0.0, // node 4 T
	  0.0, -1.0, 0.5, 1.0,			1.0, 	0.6,	0.0, // node 6
	  ]);
}

function makeFishWing() {
	fishwingVerts = new Float32Array([
	 //FISH WING
	 0.5, 0.5, 0.0, 1.0,     0.6, 0.7, 0.9,	// Node 7  PURPLE 
	 0.75, 0.25, 0.0, 1.0,   1.0, 	0.6,	0.0,	// Node 6
	 0.0, 0.0, 0.025, 1.0,     1.0, 	1.0,	0.0,	// Node 2
  
	 -0.05, 0.05, 0.025, 1.0,   1.0, 	0.6,	0.0,	// Node 1
	 0.0, 0.0, 0.025, 1.0,     0.6, 0.7, 0.9,	// Node 2
	 0.5, 0.5, 0.0, 1.0,     0.6, 0.7, 0.9,	// Node 7
  
	 // top face
	 0.5, 0.5, 0.1, 1.0,     0.6, 0.7, 0.9,	// Node 7 T -> NODE 4 PUIRPLE
	 0.75, 0.25, 0.1, 1.0,   1.0, 	0.6,	0.0,	// Node 6 T -> NODE 5
	 0.0, 0.0, 0.075, 1.0,     1.0, 	1.0,	1.0,	// Node 2 T -> NODE 3
  
	 -0.05, 0.05, 0.075, 1.0,   1.0, 	0.6,	0.0,	// Node 1 T -> NODE 0
	 0.0, 0.0, 0.075, 1.0,    1.0, 	1.0,	1.0,	// Node 2 T -> NODE 3
	 0.5, 0.5, 0.1, 1.0,     0.6, 0.7, 0.9,	// Node 7 T -> NODE 4
  
	 // side length face 1
	 0.0, 0.0, 0.075, 1.0,    1.0, 	1.0,	1.0,	// Node 2 T -> NODE 3
	 0.0, 0.0, 0.025, 1.0,     1.0, 	0.6,	0.0,	// Node 2
	 0.75, 0.25, 0.0, 1.0,   0.6, 0.7, 0.9,	// Node 6
  
	 0.75, 0.25, 0.0, 1.0,   1.0, 	0.6,	0.0,	// Node 6
	 0.0, 0.0, 0.05, 1.0,     0.6, 0.7, 0.9,	// Node 2 T -> NODE 3
	 0.75, 0.25, 0.1, 1.0,   1.0, 	0.6,	0.0,	// Node 6 T -> NODE 5
  
	 // side length face 2
	 -0.05, 0.05, 0.075, 1.0,   1.0, 	1.0,	0.0,	// Node 1 T -> NODE 0
	 -0.05, 0.05, 0.025, 1.0,   1.0, 	1.0,	1.0,	// Node 1
	 0.5, 0.5, 0.0, 1.0,     0.6, 0.7, 0.9,	// Node 7
  
	 -0.05, 0.05, 0.075, 1.0,   1.0, 	0.6,	0.0,	// Node 1 T -> NODE 0
	 0.5, 0.5, 0.0, 1.0,     0.6, 0.7, 0.9,	// Node 7
	 0.5, 0.5, 0.1, 1.0,     0.6, 0.7, 0.9,	// Node 7 T -> NODE 4
  
	 // front face wide 
	 0.5, 0.5, 0.0, 1.0,     0.6, 0.7, 0.9,	// Node 7
	 0.5, 0.5, 0.1, 1.0,     0.6, 0.7, 0.9,	// Node 7 T -> NODE 4
	 0.75, 0.25, 0.0, 1.0,   0.6, 0.7, 0.9,	// Node 6
  
	 0.5, 0.5, 0.1, 1.0,     0.6, 0.7, 0.9,	// Node 7 T -> NODE 4
	 0.75, 0.25, 0.0, 1.0,   1.0, 	1.0,	1.0,	// Node 6
	 0.75, 0.25, 0.1, 1.0,   1.0, 	0.6,	0.0,	// Node 6 T -> NODE 5
  
	 // front face small square
	 -0.05, 0.05, 0.075, 1.0,   1.0, 	1.0,	0.0,	// Node 1 T -> NODE 0
	 -0.05, 0.05, 0.025, 1.0,   1.0, 	0.6,	0.0,	// Node 1
	 0.0, 0.0, 0.025, 1.0,    1.0, 	1.0,	1.0,	// Node 2
  
	 -0.05, 0.05, 0.075, 1.0,   1.0, 	1.0,	0.0,	// Node 1 T -> NODE 0
	 0.0, 0.0, 0.025, 1.0,     0.02, 	0.9,	0.4,	// Node 2
	 0.0, 0.0, 0.075, 1.0,     1.0, 	0.6,	0.0,	// Node 2 T -> NODE 3
	]);
}

function makeAxes() {
	axesVerts = new Float32Array([
	// Drawing Axes: Draw them using gl.LINES drawing primitive;
     	// +x axis RED; +y axis GREEN; +z axis BLUE; origin: GRAY
		 0.0,  0.0,  0.0, 1.0,		0.3,  0.3,  0.3,	// X axis line (origin: gray)
		 2.6,  0.0,  0.0, 1.0,		1.0,  0.3,  0.3,	// 						 (endpoint: red)
		 
		 0.0,  0.0,  0.0, 1.0,    0.3,  0.3,  0.3,	// Y axis line (origin: white)
		 0.0,  2.6,  0.0, 1.0,		0.3,  0.3,  1.0,	//						 (endpoint: green)

		 0.0,  0.0,  0.0, 1.0,		0.3,  0.3,  0.3,	// Z axis line (origin:white)
		 0.0,  0.0,  2.6, 1.0,		0.3,  0.3,  1.0,	//
		]);
	
}

function drawAll(gl, n, modelMatrix, u_ModelMatrix, viewMatrix, u_ViewMatrix, projMatrix, u_ProjMatrix) {
//==============================================================================
  // Clear <canvas>  colors AND the depth buffer
 // gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  modelMatrix.setIdentity();    // DEFINE 'world-space' coords.

  function drawArm(x, y, z, angle1, angle2){
	modelMatrix.translate(x, y, z); // EDIT TO x, y, z

	modelMatrix.rotate(angle2, 0, 1, 0.1);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);

	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);

	modelMatrix.rotate(angle1, 0, 1, 0.1);
	modelMatrix.translate(0.5, 0.5, 1.0);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);

	modelMatrix.translate(0.0, 0.0, -2);
    modelMatrix.rotate(angle1, 0, 1, 0.1);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);

	modelMatrix.translate(0.0, 0.0, -2);
    modelMatrix.rotate(angle2, 0, 1, 0.1);	
    gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);

	modelMatrix.translate(0.0, 0.0, -2);
    modelMatrix.rotate(angle2, 0, 1, 0.1);	
    gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);

	modelMatrix.translate(0.0, 0.0, -2);
    modelMatrix.rotate(angle1, 0, 1, 0.1);	
    gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);

	modelMatrix.translate(0.0, 0.0, -2);
    modelMatrix.rotate(angle1, 0, 1, 0.1);	
    gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, tetraStart/floatsPerVertex, tetraVerts.length/floatsPerVertex);
}

function drawAnemone() {
	  //--------Draw Spinning torus
    modelMatrix.translate(0, -0.9, -0.42);
    modelMatrix.scale(0.45, 0.45, 0.85);
    modelMatrix.rotate(g_angle01, 0, 0, 1);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
    gl.drawArrays(gl.TRIANGLE_STRIP, torStart/floatsPerVertex, torVerts.length/floatsPerVertex);	

   modelMatrix.scale(0.3, 0.3, 0.15);
   count = 0;
   	for (i = 0; i < 11; i++) {
	  pushMatrix(modelMatrix);
	   drawArm(-0.4+(5.3 * Math.cos(30 + count)), -0.4 +(5.3 * Math.sin(30 + count)), -2.5, g_angle1now, g_angle0now);
	   count += 20;
	  modelMatrix = popMatrix();
   }
}

  pushMatrix(modelMatrix);  // SAVE world drawing coords.
  modelMatrix.scale(1, 1, -1);
  	drawAnemone();	
  modelMatrix = popMatrix();  // RESTORE 'world' drawing coords.

  var t = 0;
  function drawFish(){
	//---- Draw FISH BODY -
	if (t == 3) {
		quatMatrix.setFromQuat(qTot.x, qTot.y, qTot.z, qTot.w);	// Quaternion-->Matrix
		modelMatrix.concat(quatMatrix);	// apply that matrix.
	}

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.TRIANGLES, fishbodStart/floatsPerVertex, fishbodVerts.length/floatsPerVertex); 

	if (t == 3) {
		gl.drawArrays(gl.LINES, axesStart/floatsPerVertex, axesVerts.length/floatsPerVertex); // AXES
	}
	t += 1;
	

	pushMatrix(modelMatrix);
	  modelMatrix.translate(0, 0, 0.5);
	  modelMatrix.rotate(-45, 1.0, 0, 1);
	  modelMatrix.rotate(g_angle2now, 0, 1, 0.1); 
	  modelMatrix.scale(2, 2, 3);
	  gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	  gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	  gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	  gl.drawArrays(gl.TRIANGLES, fishwingStart/floatsPerVertex, fishwingVerts.length/floatsPerVertex); 
  
	  modelMatrix.translate(-0.1, 0, 0);
	  modelMatrix.rotate(180, 1.0, 0, 90);
	  modelMatrix.rotate(g_angle2now, 0, 0, 0.1);	 
	  gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	  gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	  gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	  gl.drawArrays(gl.TRIANGLES, fishwingStart/floatsPerVertex, fishwingVerts.length/floatsPerVertex); 
	modelMatrix = popMatrix();
	}

   pushMatrix(modelMatrix);
   modelMatrix.scale(0.15, 0.2, 0.15); //
   for (i = 0; i < 5; i++) {
	   if (i % 2 == 0) {
		   pushMatrix(modelMatrix);
		   modelMatrix.translate(4 * i, 0.7 + Math.pow(-1, i) * i, 0);
		   modelMatrix.rotate(40, 0, 0, 0.1);
	   }
	   else {
		   pushMatrix(modelMatrix);
		   modelMatrix.translate(0, 3 + 2 * i, 0);
	   }
	   drawFish();
	   modelMatrix = popMatrix();
	}
   modelMatrix = popMatrix();


   pushMatrix(modelMatrix);  // SAVE world drawing coords.
   //--------Draw Spinning Sphere
   modelMatrix.translate(3.0, 2.4, 0.0);
   modelMatrix.scale(0.5, 0.5, 0.5);

   gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
   gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
   gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
		   // Draw just the sphere's vertices
   gl.drawArrays(gl.TRIANGLE_STRIP, sphStart/floatsPerVertex, sphVerts.length/floatsPerVertex);	

   gl.drawArrays(gl.LINES, axesStart/floatsPerVertex, axesVerts.length/floatsPerVertex); // AXES
modelMatrix = popMatrix();  // RESTORE 'world' drawing coords.


  pushMatrix(modelMatrix);  // SAVE world drawing coords.
  	//---------Draw Ground Plane, without spinning.
  	// position it.
  	modelMatrix.translate(0.4, -0.4, 0.0);	
  	modelMatrix.scale(0.1, 0.1, 0.1);				// shrink by 10X:

  	// Drawing:
  	// Pass our current matrix to the vertex shaders:
	  gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	  gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	  gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
    // Draw just the ground-plane's vertices
    gl.drawArrays(gl.LINES, gndStart/floatsPerVertex, gndVerts.length/floatsPerVertex);

	modelMatrix.translate(-20, -20, 0);
	modelMatrix.scale(10, 10, 10);
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	gl.uniformMatrix4fv(u_ProjMatrix, false, projMatrix.elements);
	gl.uniformMatrix4fv(u_ViewMatrix, false, viewMatrix.elements);
	gl.drawArrays(gl.LINES, axesStart/floatsPerVertex, axesVerts.length/floatsPerVertex); // AXES
  modelMatrix = popMatrix();  // RESTORE 'world' drawing coords.
  //===========================================================
}
	


// Last time that this function was called:  (used for animation timing)
var g_last = Date.now();

function animate(angle) {
//==============================================================================
  // Calculate the elapsed time
  var now = Date.now();
  var elapsed = now - g_last;
  g_last = now;    
  // Update the current rotation angle (adjusted by the elapsed time)
  //  limit the angle to move smoothly between +20 and -85 degrees:
//  if(angle >  120.0 && ANGLE_STEP > 0) ANGLE_STEP = -ANGLE_STEP;
//  if(angle < -120.0 && ANGLE_STEP < 0) ANGLE_STEP = -ANGLE_STEP;

  var newAngle = angle + (g_angle01Rate * elapsed) / 1000.0;
  if(newAngle > 180.0) newAngle = newAngle - 360.0;
  if(newAngle <-180.0) newAngle = newAngle + 360.0;
  return newAngle;
 // var newAngle = angle + (ANGLE_STEP * elapsed) / 1000.0;
 // return newAngle %= 360;
}

//==================HTML Button Callbacks======================

function angleSubmit() {
	// Called when user presses 'Submit' button on our webpage

	// Read HTML edit-box contents:
	  var UsrTxt = document.getElementById('usrAngle').value;	
	// Display what we read from the edit-box: use it to fill up
	// the HTML 'div' element with id='editBoxOut':
	  document.getElementById('EditBoxOut').innerHTML ='You Typed: '+UsrTxt;
	  console.log('angleSubmit: UsrTxt:', UsrTxt); // print in console, and
	  //g_angle1now = parseFloat(UsrTxt);     // convert string to float number 
	  g_angle0now = parseFloat(UsrTxt); 
	  g_angle1now = -1 *  parseFloat(UsrTxt); 
	};
	
  
  function spinUp() {
	// Called when user presses the 'Spin >>' button on our webpage.
	  g_angle01Rate += 25; 
	}
	
  function spinDown() {
	// Called when user presses the 'Spin <<' button
	 g_angle01Rate -= 25; 
	}
  
  function runStopBase() {
	if(g_angle01Rate * g_angle01Rate > 0) {  // if nonzero rate,
	  myTmp = g_angle01Rate;  // store the current rate,
	  g_angle01Rate = 0;      // and set to zero.
	}
	else {    // but if rate is zero,
		g_angle01Rate = myTmp;  // use the stored rate.
	}
  }
  
  function runStopArms() {
	  if(g_angle0rate*g_angle0rate > 0) {  // if nonzero rate,
		myTmp1 = g_angle0rate;  // store the current rate,
		g_angle0rate = 0;      // and set to zero.
	  }
	  else {    // but if rate is zero,
		g_angle0rate = myTmp1;  // use the stored rate.
	  }
  
	  if(g_angle1rate*g_angle1rate > 0) {  // if nonzero rate,
		myTmp2 = g_angle1rate;  // store the current rate,
		g_angle1rate = 0;      // and set to zero.
	  }
	  else {    // but if rate is zero,
		g_angle1rate = myTmp2;  // use the stored rate.
	  }
	}
  
  function runStopFish() {
	  if(g_angle2rate*g_angle2rate > 0) {  // if nonzero rate,
		myTmp3 = g_angle2rate;  // store the current rate,
		g_angle2rate = 0;      // and set to zero.
	  }
	  else {    // but if rate is zero,
		g_angle2rate = myTmp3;  // use the stored rate.
	  }
	  if(g_angle3now*g_angle3now > 1) {  // if nonzero rate,
		myTmp4 = g_angle3now;  
		g_angle3now = 0;  
		g_angle3rate = 0;
		g_angle3brake = 0.0;
	  }
	  else {
		g_angle3now = myTmp4;  
		g_angle3brake = 1.0;
	  }
	}


function myMouseDown(ev, gl, canvas) {
	//==============================================================================
	// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
	  var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
	  var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
	  var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
	//  console.log('myMouseDown(pixel coords): xp,yp=\t',xp,',\t',yp);
	  
		// Convert to Canonical View Volume (CVV) coordinates too:
	  var x = (xp - canvas.width/2)  / 		// move origin to center of canvas and
							   (canvas.width/2);			// normalize canvas to -1 <= x < +1,
		var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
								 (canvas.height/2);
	//	console.log('myMouseDown(CVV coords  ):  x, y=\t',x,',\t',y);
		
		isDrag = true;											// set our mouse-dragging flag
		xMclik = x;													// record where mouse-dragging began
		yMclik = y;
	};
	
	
	function myMouseMove(ev, gl, canvas) {
	//==============================================================================
		if(isDrag==false) return;				// IGNORE all mouse-moves except 'dragging'
	
		// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
	  var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
	  var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
		var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
	//  console.log('myMouseMove(pixel coords): xp,yp=\t',xp,',\t',yp);
	  
		// Convert to Canonical View Volume (CVV) coordinates too:
	  var x = (xp - canvas.width/2)  / 		// move origin to center of canvas and
							   (canvas.width/2);			// normalize canvas to -1 <= x < +1,
		var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
								 (canvas.height/2);
	
		// find how far we dragged the mouse:
		xMdragTot += (x - xMclik);					// Accumulate change-in-mouse-position,&
		yMdragTot += (y - yMclik);
		// AND use any mouse-dragging we found to update quaternions qNew and qTot.
		dragQuat(x - xMclik, y - yMclik);
		
		xMclik = x;													// Make NEXT drag-measurement from here.
		yMclik = y;
		
		// Show it on our webpage, in the <div> element named 'MouseText':
		document.getElementById('MouseText').innerHTML=
				'Mouse Drag totals (CVV x,y coords):\t'+
				 xMdragTot.toFixed(5)+', \t'+
				 yMdragTot.toFixed(5);	
	};
	
	function myMouseUp(ev, gl, canvas) {
	//==============================================================================
	// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
	  var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
	  var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
		var yp = canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
	//  console.log('myMouseUp  (pixel coords): xp,yp=\t',xp,',\t',yp);
	  
		// Convert to Canonical View Volume (CVV) coordinates too:
	  var x = (xp - canvas.width/2)  / 		// move origin to center of canvas and
							   (canvas.width/2);			// normalize canvas to -1 <= x < +1,
		var y = (yp - canvas.height/2) /		//										 -1 <= y < +1.
								 (canvas.height/2);
	//	console.log('myMouseUp  (CVV coords  ):  x, y=\t',x,',\t',y);
		
		isDrag = false;											// CLEAR our mouse-dragging flag, and
		// accumulate any final bit of mouse-dragging we did:
		xMdragTot += (x - xMclik);
		yMdragTot += (y - yMclik);
	//	console.log('myMouseUp: xMdragTot,yMdragTot =',xMdragTot,',\t',yMdragTot);
	
		// AND use any mouse-dragging we found to update quaternions qNew and qTot;
		dragQuat(x - xMclik, y - yMclik);
	
		// Show it on our webpage, in the <div> element named 'MouseText':
		document.getElementById('MouseText').innerHTML=
				'Mouse Drag totals (CVV x,y coords):\t'+
				 xMdragTot.toFixed(5)+', \t'+
				 yMdragTot.toFixed(5);	
	};
	
	function dragQuat(xdrag, ydrag) {
	//==============================================================================
	// Called when user drags mouse by 'xdrag,ydrag' as measured in CVV coords.
	// We find a rotation axis perpendicular to the drag direction, and convert the 
	// drag distance to an angular rotation amount, and use both to set the value of 
	// the quaternion qNew.  We then combine this new rotation with the current 
	// rotation stored in quaternion 'qTot' by quaternion multiply.  
	// Note the 'draw()' function converts this current 'qTot' quaternion to a  
	// rotation matrix for drawing. 
		var res = 5;
		var qTmp = new Quaternion(0,0,0,1);
		
		var dist = Math.sqrt(xdrag*xdrag + ydrag*ydrag);
		qNew.setFromAxisAngle(-ydrag + 0.0001, xdrag + 0.0001, 0.0, dist*150.0);
								
		qTmp.multiply(qNew,qTot);			// apply new rotation to current rotation. 
		// Matrix4.prototype.setFromQuat().
		//qTmp.normalize();						// normalize to ensure we stay at length==1.0.
		qTot.copy(qTmp);
		// show the new quaternion qTot on our webpage in the <div> element 'QuatValue'
	};

	function keydown(ev, gl) {
		if(ev.keyCode == 39) { // The right arrow key was pressed, we want to strafe right
			eyeX += 0.1 * Math.sin(thetarad);
			eyeY -= 0.1 * Math.cos(thetarad);
		} 
		else 
		if (ev.keyCode == 37) { // The left arrow key was pressed
			eyeX -= 0.1 * Math.sin(thetarad);
			eyeY += 0.1 * Math.cos(thetarad);
		} 
		else 
		if (ev.keyCode == 38) { // Up arrow key
			eyeX += 0.1 * Math.cos(thetarad);
			eyeY += 0.1 * Math.sin(thetarad);
			eyeZ += 0.1 * deltaz;
		} 

		else
		if (ev.keyCode == 40) { // Down arrow key
			eyeX -= 0.1 * Math.cos(thetarad);
			eyeY -= 0.1 * Math.sin(thetarad);
			eyeZ -= 0.1 * deltaz;
		} 

		else
		if (ev.keyCode == 87) { // W key
			deltaz += 0.05;
		} 
		else
		if (ev.keyCode == 65) { // A key
			thetarad += Math.PI/100;
		} 
		else
		if (ev.keyCode == 83) { // S key
			deltaz -= 0.05;
		} 
		else
		if (ev.keyCode == 68) { // D key
			thetarad -= Math.PI/100;
		}
		
		else { return; }

		drawAll(gl, n, modelMatrix, u_ModelMatrix, viewMatrix, u_ViewMatrix, projMatrix, u_ProjMatrix); 
}


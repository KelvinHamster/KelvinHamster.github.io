

{
    
    let n2container = document.getElementById('numericalcanvas2');
    let n2 = n2container.getContext('2d');

    
    let n2container3D = document.getElementById('numericalcanvas3');
    let n2_3D = n2container3D.getContext('2d');


    let N = 10;
    let M = 10;
    let h = 1;
    let R1 = 2;
    let R2 = 2;
    let dt = 0.02;

    let t = 0;


    let dot_size = 0.1/N;
    let stroke_size_2d = 5;
    let stroke_size_3d = 5;


    let x0_1 = 0;
    let x0_2 = 0;
    let a_1 = 1;
    let a_2 = 1;
    let X = new Float64Array((N+1)*M*6);
    //X[m,n] = [x,y,z,x',y',z']
    
    let X_ = new Float64Array((N+1)*M*6);//temporary for swapping

    let beta = 1;
    let alpha = 1;

    let IC_dev = 0;

    var RUN_N2 = false;


    let init_X = function(){
        //solve the catenoid
        var cats = find_catenoids(h/R1,R2/R1);
        a_1 = cats[0]*R1;
        a_2 = cats[2]*R1;
        x0_1 = cats[1]*a_1;
        x0_2 = cats[3]*a_2;
        num_steps = 0;

        
        X = new Float64Array((N+1)*M*6);
        X_ = new Float64Array((N+1)*M*6);
        
        for(let m = 0; m < M; m++){
            let index = 6*(m*(N+1));
            let angle = (2*Math.PI/M)*m;
            X[index] = 0;
            X[index+1] = R1 * Math.cos(angle);
            X[index+2] = R1 * Math.sin(angle);
            X[index+3] = 0; X[index+4] = 0; X[index+5] = 0;
            X_[index] = 0;
            X_[index+1] = R1 * Math.cos(angle);
            X_[index+2] = R1 * Math.sin(angle);
            X_[index+3] = 0; X_[index+4] = 0; X_[index+5] = 0;
            index = 6*(m*(N+1) + N);
            X[index] = h;
            X[index+1] = R2 * Math.cos(angle);
            X[index+2] = R2 * Math.sin(angle);
            X[index+3] = 0; X[index+4] = 0; X[index+5] = 0;
            X_[index] = h;
            X_[index+1] = R2 * Math.cos(angle);
            X_[index+2] = R2 * Math.sin(angle);
            X_[index+3] = 0; X_[index+4] = 0; X_[index+5] = 0;
        }
        for(let k = 1; k < N; k++){
            let t = k/N;
            let r = 0;
            switch(init_cond){
                case 0:
                    r = (1-t)*R1 + t*R2;
                    break;
                case 1:
                    r = a_1*Math.cosh((t*h-x0_1)/a_1);
                    break;
                case 2:
                    r = a_2*Math.cosh((t*h-x0_2)/a_2);
            }
            r *= (1+IC_dev);
            for(let m = 0; m < M; m++){
                let index = 6*(m*(N+1) + k);
                let angle = (2*Math.PI/M)*m;
                X[index] = t*h;
                X[index+1] = r * Math.cos(angle);
                X[index+2] = r * Math.sin(angle);
                X[index+3] = 0; X[index+4] = 0; X[index+5] = 0;
            }
        }

        
        t = 0;
        update_feedback_params();
    }

    let update_feedback_params = function(){
        document.getElementById("N2_STEP_COUNT").innerHTML = t.toPrecision(3);
        let A = Area(), T = KE();
        document.getElementById("N2_A_VALUE").innerHTML = A.toPrecision(7);
        document.getElementById("N2_KINETIC").innerHTML = T.toExponential(3);
        document.getElementById("N2_ENERGY").innerHTML = (T + alpha*A).toPrecision(7);

    }

    let space_transform = function(x,y,z){
        //3D->2D using (x,r) coordinates
        y = Math.sqrt(y*y + z*z);

        //convert to canvas space
        let maxr = Math.max(R1,R2);
        let margin = 20;
        let width = n2container.width - margin*2
        let height = n2container.height - margin*2;
        //scale factor:
        let scale_factor = Math.min(width/h, height/maxr);
        x = margin + x*scale_factor;
        y = margin + height - y*scale_factor;
        return Array(x,y);
    }

    let space_transform_3D = function(x,y,z){
        //3D->2D using perspective transform
        /*
          simple way is to stand at x=h/2, y=0, z= -D
          looking in direction +z, so we project onto z=-D + 1

          however, for something more visual, do some rotations. the
          transformations will thus be

          rotate about x-axis by theta -> translate so (h/2,0,0) is origin
          -> rotate about new y-axis by theta -> perspective project

          Some information on the perspective projection is found here
          https://en.wikipedia.org/wiki/Transformation_matrix#Perspective_projection
        */
        let theta = 0.3;
        let costheta = Math.cos(theta),sintheta = Math.sin(theta);
        [y,z] = [y*costheta - z*sintheta, y*sintheta + z*costheta];
        x -= h/2;
        [x,z] = [x*costheta - z*sintheta, x*sintheta + z*costheta];

        let maxr = Math.max(R1,R2);
        let D = -2*maxr;
        z += D;

        x /= z;
        y/= z;
        //z /= y;

        //convert to canvas space
        let margin = 20;
        let width = n2container3D.width - margin*2
        let height = n2container3D.height - margin*2;
        //scale factor:
        let scale_factor = Math.min(width/h, height/maxr);
        x = margin + width/2 + x*scale_factor;
        y = margin + height/2 - y*scale_factor;
        return Array(x,y);
    }

    /**
     * Calculates the area of the triangle that comes from the
     * convex hull of the origin, v1, and v2
     */
    let triangle_calc = function(v1x,v1y,v1z,v2x,v2y,v2z){
        //cross product
        let x = v1y * v2z - v1z * v2y;
        let y = v1z * v2x - v1x * v2z;
        let z = v1x * v2y - v1y * v2x;
        return Math.sqrt(x*x + y*y + z*z)/2;
    }

    /**
     * Calculates the partial derivative of A(v1,v2) with respect
     * to r, where v1 = r - r_1, v2 = r - r_2,
     * and add them to the array 'total'.
     */
    let partial_traingle_calc = function(total,v1x,v1y,v1z,v2x,v2y,v2z){
        //cross product (v1 x v2)
        let v1v2x = v1y * v2z - v1z * v2y;
        let v1v2y = v1z * v2x - v1x * v2z;
        let v1v2z = v1x * v2y - v1y * v2x;

        let denom = 2*Math.sqrt(v1v2x*v1v2x + v1v2y*v1v2y + v1v2z*v1v2z)

        //cross product (e_0 x v) = (0,-z,y)
        total[0] += ((v2y-v1y)*v1v2z - (v2z-v1z)*v1v2y)/denom;
        //cross product (e_1 x v) = (z,0,-x)
        total[1] += ((v2z-v1z)*v1v2x - (v2x-v1x)*v1v2z)/denom;
        //cross product (e_2 x v) = (-y,x,0)
        total[2] += ((v2x-v1x)*v1v2y - (v2y-v1y)*v1v2x)/denom;
    }

    /**
     * Calculates the area function for a given X
     */
     let Area = function(){
        /*
          we are using a rectangular lattice, so area is not very well
          defined. Approximate using two triangles.
        */
        let A  = 0;
        let bigmod = (N+1)*M*6;
        let m_shift = (N+1)*6;
        
        for(let m = 0; m < M; m++){
            let index0 = m*(N+1);
            for(let n = 0; n < N; n++){
                let index = 6*(index0 + n);
                //index + 6 => X_{n+1,m}
                //index + 6*(N+1) => X_{n,m+1}
                let minc_index = (index+m_shift) % bigmod//index of X_{n,m+1}
                //this will be the same for both triangles
                //v2 = X_{n+1,m+1} - X_{n_m}
                let v2x = X[minc_index + 6] - X[index];
                let v2y = X[minc_index + 7] - X[index+1];
                let v2z = X[minc_index + 8] - X[index+2];
                A += triangle_calc(
                    X[index+6] - X[index],//v1 = X_{n+1,m} - X_{n_m}
                    X[index+7] - X[index+1],
                    X[index+8] - X[index+2],
                    v2x,v2y,v2z
                ) + triangle_calc(
                    X[minc_index] - X[index],//v1 = X_{n,m+1} - X_{n_m}
                    X[minc_index+1] - X[index+1],
                    X[minc_index+2] - X[index+2],
                    v2x,v2y,v2z
                );
            }
        }
        return A;
    }

    let timestep = function(){
        let bigmod = (N+1)*M*6;
        let m_shift = (N+1)*6;

        let deriv_sum = new Float64Array(3);
        for(let m = 0; m < M; m++){
            let index0 = m*(N+1);
            for(let n = 1; n < N; n++){
                let index = 6*(index0 + n);
                let minc = (index+m_shift) % bigmod//index of X_{n,m+1}
                let msub = (index+bigmod-m_shift) % bigmod//index of X_{n,m-1}
                //euler's
                //position
                X_[index] = X[index] + X[index+3]*dt;
                X_[index+1]=X[index+1]+X[index+4]*dt;
                X_[index+2]=X[index+2]+X[index+5]*dt;

                //velocity
                /*
                  six triangles to worry about. For each, the other vertices
                  correspond to x_ ...
                  (n-1,m-1), (n-1,m)
                  (n-1,m-1), (n,m-1)
                  (n+1,m+1), (n+1,m)
                  (n+1,m+1), (n,m+1)
                  (n+1,m)  , (n,m-1)
                  (n-1,m)  , (n,m+1)
                  to solve for the derivatives of their areas, utilize triangle
                  symmetries
                */
                deriv_sum[0] = 0; deriv_sum[1] = 0; deriv_sum[2] = 0;
                let x = X[index], y = X[index+1], z=X[index+2];
                partial_traingle_calc(deriv_sum,
                    x-X[msub -6], y-X[msub -5], z-X[msub -4],
                    x-X[index-6], y-X[index-5], z-X[index-4]);
                partial_traingle_calc(deriv_sum,
                    x-X[msub-6], y-X[msub-5], z-X[msub-4],
                    x-X[msub  ], y-X[msub+1], z-X[msub+2]);
                partial_traingle_calc(deriv_sum,
                    x-X[minc +6], y-X[minc +7], z-X[minc +8],
                    x-X[index+6], y-X[index+7], z-X[index+8]);
                partial_traingle_calc(deriv_sum,
                    x-X[minc+6], y-X[minc+7], z-X[minc+8],
                    x-X[minc  ], y-X[minc+1], z-X[minc+2]);
                partial_traingle_calc(deriv_sum,
                    x-X[index+6], y-X[index+7], z-X[index+8],
                    x-X[msub   ], y-X[msub +1], z-X[msub +2]);
                partial_traingle_calc(deriv_sum,
                    x-X[index-6], y-X[index-5], z-X[index-4],
                    x-X[minc   ], y-X[minc +1], z-X[minc +2]);
                X_[index+3]=(1-2*beta*dt)*X[index+3]-dt*alpha*deriv_sum[0];
                X_[index+4]=(1-2*beta*dt)*X[index+4]-dt*alpha*deriv_sum[1];
                X_[index+5]=(1-2*beta*dt)*X[index+5]-dt*alpha*deriv_sum[2];
            }
        }
        //swap X and X_
        [X, X_] = [X_, X];
        t += dt;
    }

    let KE = function(){
        let total = 0;
        let m_shift = (N+1)*6;
        for(let n = 1; n < N; n++){
            let offset = n*6 + 3;
            for(let m = 0; m < M; m++){
                total += X[offset]*X[offset] +
                    X[offset+1]*X[offset+1] +
                    X[offset+2]*X[offset+2];
                offset += m_shift;
            }
        }
        return total/2;
    }

    var stepN2 = function(){
        timestep();
        updateN2();
        update_feedback_params();
    }

    var updateN2 = function(){
        n2.clearRect(0, 0, n2container.width, n2container.height);
        n2_3D.clearRect(0, 0, n2container3D.width, n2container3D.height);
        //draw catenoids
        if(!isNaN(a_1)){
            n2.lineWidth = stroke_size_2d/2;
            n2.strokeStyle = 'rgba(0,0,180,0.5)';
            n2.beginPath();
            n2.moveTo(...space_transform(0,R1,0));
            let plotstep = 0.2 * h/N;
            for(let t = 0; t < h; t+= plotstep){
                n2.lineTo(...space_transform(t,a_1*Math.cosh((t-x0_1)/a_1),0));
            }
            n2.lineTo(...space_transform(h,R2,0));
            n2.stroke();
            n2.strokeStyle = 'rgba(0,180,0,0.5)';
            n2.beginPath();
            n2.moveTo(...space_transform(0,R1,0));
            for(let t = 0; t < h; t+= plotstep){
                n2.lineTo(...space_transform(t,a_2*Math.cosh((t-x0_2)/a_2),0));
            }
            n2.lineTo(...space_transform(h,R2,0));
            n2.stroke();
        }
        

        n2.lineWidth = stroke_size_2d;
        n2_3D.lineWidth = stroke_size_3d;
        n2.fillStyle = 'rgb(230,150,20)';
        n2.strokeStyle = 'rgb(0,0,0)';
        n2_3D.strokeStyle = 'rgb(0,0,0)';
        n2.beginPath();
        n2_3D.beginPath();
        //boundary
        n2.moveTo(...space_transform(0,R1,0));
        n2.lineTo(...space_transform(0,0,0));
        n2.lineTo(...space_transform(h,0,0));
        n2.lineTo(...space_transform(h,R2,0));
        
        n2_3D.moveTo(...space_transform_3D(0,R1,0));
        n2_3D.lineTo(...space_transform_3D(0,-R1,0));
        n2_3D.moveTo(...space_transform_3D(h,R2,0));
        n2_3D.lineTo(...space_transform_3D(h,-R2,0));
        n2_3D.moveTo(...space_transform_3D(0,0,R1));
        n2_3D.lineTo(...space_transform_3D(0,0,-R1));
        n2_3D.moveTo(...space_transform_3D(h,0,R2));
        n2_3D.lineTo(...space_transform_3D(h,0,-R2));
        n2_3D.moveTo(...space_transform_3D(0,0,0));
        n2_3D.lineTo(...space_transform_3D(h,0,0));
        //ticks
        let tick = 0.5;
        let tickwidth = 0.05;
        for(let y = tick; y< R1; y+= tick){
            n2.moveTo(...space_transform(0,y,0));
            n2.lineTo(...space_transform(tickwidth,y,0));
        }
        for(let y = tick; y< R2; y+= tick){
            n2.moveTo(...space_transform(h,y,0));
            n2.lineTo(...space_transform(h-tickwidth,y,0));
        }
        for(let x = tick; x< h; x += tick){
            n2.moveTo(...space_transform(x,0,0));
            n2.lineTo(...space_transform(x,tickwidth,0));
        }
        n2.stroke();
        n2_3D.stroke();
        
        n2.strokeStyle = 'rgb(0,0,200)';
        n2_3D.strokeStyle = 'rgb(0,0,200)';
        n2.beginPath();
        n2_3D.beginPath();
        for(let m = 0; m < M; m++){
            let offset = m*(N+1)*6;
            n2.moveTo(...space_transform(0,R1,0));
            n2_3D.moveTo(...space_transform_3D(X[offset],X[offset+1],X[offset+2]));
            for(let k = 1; k <= N; k++){
                let offset2 = offset + k*6;
                n2.lineTo(...space_transform(X[offset2],X[offset2+1],X[offset2+2]));
                n2_3D.lineTo(...space_transform_3D(X[offset2],X[offset2+1],X[offset2+2]));
            }
        }
        n2.stroke();
        let m_shift = (N+1)*6;
        for(let n = 1; n < N; n++){
            let offset = n*6;
            n2_3D.moveTo(...space_transform_3D(X[offset],X[offset+1],X[offset+2]));
            for(let m = 0; m < M; m++){
                offset += m_shift;
                n2_3D.lineTo(...space_transform_3D(X[offset],X[offset+1],X[offset+2]));
            }
            offset = n*6;
            n2_3D.lineTo(...space_transform_3D(X[offset],X[offset+1],X[offset+2]));
        }
        n2_3D.stroke()
        

        //no vertices on 3d because that will look weird
        n2.beginPath();
        for(let m = 0; m < M; m++){
            let offset = m*(N+1)*6;
            for(let k = 1; k < N; k++){
                let offset2 = offset + k*6;
                n2.moveTo(...space_transform(X[offset2],X[offset2+1],X[offset2+2]));
                n2.arc(...space_transform(X[offset2],X[offset2+1],X[offset2+2]),
                    n2container.width*dot_size, 0, 2 * Math.PI);
            }
        }
        n2.fill();
    }

    function N2_update_R1(val){
        R1 = parseFloat(val);
        document.getElementById("N2_SLIDER_R1").innerHTML = `R1 = ${val}`
        init_X();
        updateN2();
    };
    function N2_update_R2(val){
        R2 = parseFloat(val);
        document.getElementById("N2_SLIDER_R2").innerHTML = `R2 = ${val}`
        init_X();
        updateN2();
    };
    function N2_update_h(val){
        h = parseFloat(val);
        document.getElementById("N2_SLIDER_h").innerHTML = `h = ${val}`
        init_X();
        updateN2();
    };
    function N2_update_N(val){
        N = Math.round(Math.pow(10,parseFloat(val)));
        document.getElementById("N2_SLIDER_N").innerHTML = `N = ${N}`
        dot_size = 0.1/N;
        init_X();
        updateN2();
    };
    function N2_update_M(val){
        M = Math.round(Math.pow(10,parseFloat(val)));
        document.getElementById("N2_SLIDER_M").innerHTML = `M = ${M}`
        init_X();
        updateN2();
    };
    function N2_update_IC_dev(val){
        IC_dev = parseFloat(val);
        if(IC_dev < 0){
            document.getElementById("N2_SLIDER_IC_DEV").innerHTML =
                    `Initial Condition offset = -${(-100*IC_dev).toPrecision(2)}%`
        }else{
            document.getElementById("N2_SLIDER_IC_DEV").innerHTML =
                    `Initial Condition offset = +${(100*IC_dev).toPrecision(2)}%`
        }
        init_X();
        updateN2();
    }
    function N2_toggle_run(state){
        RUN_N2 = state;
    }
    function N2_set_initial_condition(val){
        init_cond = parseInt(val);
        init_X();
        updateN2();
    }
    function N2_update_beta(val){
        beta = parseFloat(val);
        document.getElementById("N2_SLIDER_beta").innerHTML = `beta = ${val}`
    };
    function N2_update_dt(val){
        dt = parseFloat(val);
        document.getElementById("N2_SLIDER_dt").innerHTML = `time step = ${val}`
    };
    

    let resize_canvas = function(){
        n2container.width  = window.innerWidth;
        n2container.height = n2container.width*0.6;
        n2container3D.width  = window.innerWidth;
        n2container3D.height = n2container3D.width*0.6;
        updateN2();
    }
    window.addEventListener('resize', resize_canvas, false);
    init_X();
    resize_canvas();

    function N2_estimate_a(){
        //calculate s0; use one sided difference
        let s0 = (Math.sqrt(X[7]*X[7] + X[8]*X[8]) - R1)/(X[6]);
        //B/A = -sinh^{-1}(s0),   A = R1/cosh(-B/A)
        let BA = -Math.asinh(s0);
        let A = R1/Math.cosh(-BA);
        return A;
    }
}
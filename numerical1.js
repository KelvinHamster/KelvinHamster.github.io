{
    let n1container = document.getElementById('numericalcanvas1');
    let n1 = n1container.getContext('2d');
    
    let N = 10;
    let r = new Float64Array(N+1);
    let h = 1;
    let R1 = 2;
    let R2 = 2;
    let dot_size = 0.1/N;
    let stroke_size = 5;

    let x0_1 = 0;
    let x0_2 = 0;
    let a_1 = 1;
    let a_2 = 1;
    let num_steps = 0;
    let IC_dev = 0;

    var RUN_N1 = false;
    var init_cond = 0;
    
    let init_r = function(){
        //solve the catenoid
        var cats = find_catenoids(h/R1,R2/R1);
        a_1 = cats[0]*R1;
        a_2 = cats[2]*R1;
        x0_1 = cats[1]*a_1;
        x0_2 = cats[3]*a_2;
        num_steps = 0;

        r[0] = R1; r[N] = R2;
        for(let k = 1; k < N; k++){
            let t = k/N;
            switch(init_cond){
                case 0:
                    r[k] = (1-t)*R1 + t*R2;
                    break;
                case 1:
                    r[k] = a_1*Math.cosh((t*h-x0_1)/a_1);
                    break;
                case 2:
                    r[k] = a_2*Math.cosh((t*h-x0_2)/a_2);
            }
            r[k] *= (1+IC_dev);
        }
        
        let grad = gradient(r);
        document.getElementById("N1_GRAD_NORM").innerHTML =
            Math.sqrt(vec_dot(grad,grad)).toExponential(4);
        document.getElementById("N1_STEP_COUNT").innerHTML = num_steps;
        document.getElementById("N1_A_VALUE").innerHTML = Area(r);
        second_deriv_analysis();
    }
    
    let space_transform = function(x,y){
        //convert to canvas space
        let maxr = Math.max(R1,R2);
        let margin = 20;
        let width = n1container.width - margin*2
        let height = n1container.height - margin*2;
        //scale factor:
        let scale_factor = Math.min(width/h, height/maxr);
        x = margin + x*scale_factor;
        y = margin + height - y*scale_factor;
        return Array(x,y);
    }

    /**
     * Calculates the area function for a given r.
     */
    let Area = function(r){
        let A  = 0;
        let scale = h*h/(N*N);
        for(let k = 0; k < N; k++){
            A += Math.sqrt(scale + (r[k+1] - r[k])*(r[k+1] - r[k])) *
                    (r[k+1] + r[k])/2;
        }
        return A * (2*Math.PI);
    }

    /**
     * Calculates the gradient of the area function for a given r.
     */
    let gradient = function(r){
        let scale = h*h/(N*N);

        var grad = new Float64Array(N+1);
        grad[0] = 0; grad[N] = 0; //ends are fixed
        let coef = 2*Math.PI*h/N;
        for(let k = 1; k < N; k++){
            grad[k] = (r[k]*r[k] - r[k]*r[k+1] + scale/2) /
                Math.sqrt(scale + (r[k+1] - r[k])*(r[k+1] - r[k]));

            grad[k] += (r[k]*r[k] - r[k]*r[k-1] + scale/2) /
                Math.sqrt(scale + (r[k-1] - r[k])*(r[k-1] - r[k]));

            grad[k] *= coef;
        }
        return grad;
    }
    /**
     * Adds the two arrays element-wise,
     * creating a new array without modifying
     * the old one
     */
    let vec_sum = function(a,b){
        var sum = new Float64Array(N+1);
        for(let k = 0; k <= N; k++){
            sum[k] = a[k] + b[k];
        }
        return sum;
    }
    /**
     * Adds the two arrays element-wise,
     * modifying the vector a
     */
    let vec_add = function(a,b){
        for(let k = 0; k <= N; k++){
            a[k] += b[k];
        }
    }
    /**
     * Multiply the vector v by a scalar,
     * creating a new array without modifying
     * the old one
     */
     let vec_prod = function(a,v){
        var prod = new Float64Array(N+1);
        for(let k = 0; k <= N; k++){
            prod[k] = a * v[k];
        }
        return prod;
    }
    /**
     * Multiply the vector v by a scalar,
     * modifying the vector
     */
    let vec_mult = function(a,v){
        for(let k = 0; k <= N; k++){
            v[k] *= a;
        }
    }
    /**
     * Calculates the dot product between a and b
     */
    let vec_dot = function(a,b){
        var total = 0;
        for(let k = 0; k <= N; k++){
            total += a[k] * b[k];
        }
        return total;
    }
    let gradient_step = function(){
        var grad = gradient(r);
        var eta = 1/2;
        var gamma = 1/2;
        var alpha0 = R1/Math.sqrt(vec_dot(grad,grad));
        var Ar = Area(r);
        var p = vec_prod(-alpha0,grad);//this will encompass alpha_0 gamma^i p
        while(Area(vec_sum(r,p)) - Ar > eta * vec_dot(grad,p)){
            vec_mult(gamma,p);
        }
        //step
        vec_add(r,p);

        //force r > 0
        for(let k = 0; k <= N; k++){
            r[k] = Math.max(r[k],0);
        }
    }

    var stepN1 = function(){
        gradient_step();
        updateN1();
        num_steps++;
        let grad = gradient(r);
        document.getElementById("N1_GRAD_NORM").innerHTML =
            Math.sqrt(vec_dot(grad,grad)).toExponential(4);
        document.getElementById("N1_STEP_COUNT").innerHTML = num_steps;
        document.getElementById("N1_A_VALUE").innerHTML = Area(r);
        second_deriv_analysis();
    }
    
    var updateN1 = function(){
        n1.clearRect(0, 0, n1container.width, n1container.height);
        //draw catenoids
        if(!isNaN(a_1)){
            n1.lineWidth = stroke_size/2;
            n1.strokeStyle = 'rgba(0,0,180,0.5)';
            n1.beginPath();
            n1.moveTo(...space_transform(0,R1));
            let plotstep = 0.2 * h/N;
            for(let t = 0; t < h; t+= plotstep){
                n1.lineTo(...space_transform(t,a_1*Math.cosh((t-x0_1)/a_1)));
            }
            n1.lineTo(...space_transform(h,R2));
            n1.stroke();
            n1.strokeStyle = 'rgba(0,180,0,0.5)';
            n1.beginPath();
            n1.moveTo(...space_transform(0,R1));
            for(let t = 0; t < h; t+= plotstep){
                n1.lineTo(...space_transform(t,a_2*Math.cosh((t-x0_2)/a_2)));
            }
            n1.lineTo(...space_transform(h,R2));
            n1.stroke();
        }
        

        n1.lineWidth = stroke_size;
        n1.fillStyle = 'rgb(230,150,20)';
        n1.strokeStyle = 'rgb(0,0,0)';
        n1.beginPath();
        //boundary
        n1.moveTo(...space_transform(0,R1));
        n1.lineTo(...space_transform(0,0));
        n1.lineTo(...space_transform(h,0));
        n1.lineTo(...space_transform(h,R2));
        //ticks
        let tick = 0.5;
        let tickwidth = 0.05;
        for(let y = tick; y< R1; y+= tick){
            n1.moveTo(...space_transform(0,y));
            n1.lineTo(...space_transform(tickwidth,y));
        }
        for(let y = tick; y< R2; y+= tick){
            n1.moveTo(...space_transform(h,y));
            n1.lineTo(...space_transform(h-tickwidth,y));
        }
        for(let x = tick; x< h; x += tick){
            n1.moveTo(...space_transform(x,0));
            n1.lineTo(...space_transform(x,tickwidth));
        }
        n1.stroke();
        
        n1.strokeStyle = 'rgb(0,0,200)';
        n1.beginPath();
        n1.moveTo(...space_transform(0,R1));
        for(let k = 1; k < N; k++){
            let t = k*h/N;
            n1.lineTo(...space_transform(t,r[k]));
        }
        n1.lineTo(...space_transform(h,R2));
        n1.stroke()
    
        n1.beginPath();
        for(let k = 1; k < N; k++){
            let t = k*h/N;
            n1.moveTo(...space_transform(t,r[k]))
            n1.arc(...space_transform(t,r[k]),
                n1container.width*dot_size, 0, 2 * Math.PI);
        }
        n1.fill();
    }

    function N1_update_R1(val){
        R1 = parseFloat(val);
        document.getElementById("N1_SLIDER_R1").innerHTML = `R1 = ${val}`
        init_r();
        updateN1();
    };
    function N1_update_R2(val){
        R2 = parseFloat(val);
        document.getElementById("N1_SLIDER_R2").innerHTML = `R2 = ${val}`
        init_r();
        updateN1();
    };
    function N1_update_h(val){
        h = parseFloat(val);
        document.getElementById("N1_SLIDER_h").innerHTML = `h = ${val}`
        init_r();
        updateN1();
    };
    function N1_update_N(val){
        N = Math.round(Math.pow(10,parseFloat(val)));
        document.getElementById("N1_SLIDER_N").innerHTML = `N = ${N}`
        r = new Float64Array(N+1);
        dot_size = 0.1/N;
        init_r();
        updateN1();
    };
    function N1_update_IC_dev(val){
        IC_dev = parseFloat(val);
        if(IC_dev < 0){
            document.getElementById("N1_SLIDER_IC_DEV").innerHTML =
                    `Initial Condition offset = -${(-100*IC_dev).toPrecision(2)}%`
        }else{
            document.getElementById("N1_SLIDER_IC_DEV").innerHTML =
                    `Initial Condition offset = +${(100*IC_dev).toPrecision(2)}%`
        }
        init_r();
        updateN1();
    }
    function N1_toggle_run(state){
        RUN_N1 = state;
    }
    function N1_set_initial_condition(val){
        init_cond = parseInt(val);
        init_r();
        updateN1();
    }

    let second_partial = function(i,j){
        //multivariate finite difference
        let EPS = 0.00001
        rpp = vec_prod(1,r); rpp[i] += EPS; rpp[j] += EPS;
        rpm = vec_prod(1,r); rpm[i] += EPS; rpm[j] -= EPS;
        rmp = vec_prod(1,r); rmp[i] -= EPS; rmp[j] += EPS;
        rmm = vec_prod(1,r); rmm[i] -= EPS; rmm[j] -= EPS;
        return (Area(rpp) - Area(rpm) - Area(rmp) + Area(rmm))/(4*EPS*EPS)
    }
    let second_deriv_analysis = function(){
        let maxD = minD = second_partial(1,1);//min/max second derivatives
        let maxI = minI = 1;//indices of above
        for(let k = 2; k < N ; k++){
            let deriv = second_partial(k,k);
            if(deriv > maxD){
                maxD = deriv; maxI = k;
            }else if(deriv < minD){
                minD = deriv; minI = k;
            }
        }
        let v = new Float64Array(N+1);for(k=1;k<N;k++){v[k] = 1};
        document.getElementById("N1_HESS_DIAG_MAX").innerHTML =
            `Node ${maxI}, Value: ${maxD.toExponential(3)}`;
        document.getElementById("N1_HESS_DIAG_MIN").innerHTML =
            `Node ${minI}, Value: ${minD.toExponential(3)}`;
        document.getElementById("N1_HESS_D2_same").innerHTML =
            `Value: ${second_directional_derivative(v)}`;
    }
    let second_directional_derivative = function(v){
        let total = 0;
        for(i = 1; i < N; i++){
            let row = 0;
            for(j = 1; j < N; j++){
                row += second_partial(i,j)*v[j];
            }
            total += row * v[i];
        }
        return total;
    }

    let resize_canvas = function(){
        n1container.width  = window.innerWidth;
        n1container.height = n1container.width*0.6;
        updateN1();
    }
    window.addEventListener('resize', resize_canvas, false);
    init_r();
    resize_canvas();


    /* interfacing functions for data collection */
    function N1_estimate_a(){
        //calculate s0; use one sided difference
        let s0 = (r[1]-r[0])/(h/N);
        //B/A = -sinh^{-1}(s0),   A = R1/cosh(-B/A)
        let BA = -Math.asinh(s0);
        let A = R1/Math.cosh(-BA);
        return A;
    }
}

/* NOT SPECIFIC TO NUMERICAL1 */


function find_catenoids(h,alpha,tol = 0.0001){
    function Fp(a){
        let temp = Math.sqrt(1 - a*a);
        let temp2 = 1/(a*a);
        return (a/temp - h*temp2)*Math.sinh(h/a) + temp2*temp*h*Math.cosh(h/a);
    }
    function F(a){
        return Math.cosh(h/a) - Math.sinh(h/a)*Math.sqrt(1 - a*a) - alpha;
    }

    //Find a_ that has F(a_) < 0 using regular gradient descent
    let a_ = 0.5;
    let grad = 0, Fa = 0, eta = 0, gamma = 0, alpha0 = 0, p = 0;
    let max_iters = 100;
    while(F(a_) > 0 && max_iters > 0){
        grad = Fp(a_);
        if(Math.abs(grad) < 0.000001){
            //safe to assume no solution exists
            return [NaN,NaN,NaN,NaN];
        }
        Fa = F(a_);
        gamma = 0.5;
        eta = 0.5;

        alpha0 = 0.7*Math.min(a_,1-a_);
        p = -alpha0*Math.sign(grad);
        while(F(a_ + p) - Fa > eta*grad*p){
            p *= gamma;
        }
        a_ += p;
        max_iters --;
    }
    if(max_iters == 0){
        //safe to assume no solution exists
        return [NaN,NaN,NaN,NaN];
    }

    //negative a_, so modified bisection method for first
    let bis_a = 0, bis_b = a_, bis_c = a_/2;
    //bis_a will evaluate to infty or > 0, bis_b will evaluate to < 0
    let Fc = F(bis_c);//midpoint value
    while(bis_b - bis_a > tol){
        if(Fc < 0){
            //zero is between a and c
            bis_b = bis_c;
            bis_c = (bis_a + bis_b)/2;
            Fc = F(bis_c);
        }else{
            //zero is between c and b
            bis_a = bis_c;
            bis_c = (bis_a + bis_b)/2;
            Fc = F(bis_c);
        }
    }
    // result array by (a1,x0_1/R_1,a2,x0_2/R_1)
    var RESULT = [bis_c,Math.acosh(1/bis_c),0,0];
    //other value?
    if(Math.cosh(h) > alpha){
        //positive h: do bisection method again
        bis_a = a_; bis_b = 1; bis_c = (1+a_)/2;
        //F(bis_a) < 0, F(bis_b) > 0
        Fc = F(bis_c);//midpoint value
        while(bis_b - bis_a > tol){
            if(Fc > 0){
                //zero is between a and c
                bis_b = bis_c;
                bis_c = (bis_a + bis_b)/2;
                Fc = F(bis_c);
            }else{
                //zero is between c and b
                bis_a = bis_c;
                bis_c = (bis_a + bis_b)/2;
                Fc = F(bis_c);
            }
        }
        RESULT[2] = bis_c;
        RESULT[3] = Math.acosh(1/bis_c);
    }else{
        //negative h: bisect between 0 and 1.
        h = -h;
        bis_a = 0; bis_b = 1; bis_c = 0.5;
        //F(bis_a) > 0, F(bis_b) < 0
        Fc = F(bis_c);//midpoint value
        while(bis_b - bis_a > tol){
            if(Fc < 0){
                //zero is between a and c
                bis_b = bis_c;
                bis_c = (bis_a + bis_b)/2;
                Fc = F(bis_c);
            }else{
                //zero is between c and b
                bis_a = bis_c;
                bis_c = (bis_a + bis_b)/2;
                Fc = F(bis_c);
            }
        }
        RESULT[2] = bis_c;
        RESULT[3] = -Math.acosh(1/bis_c);
    }
    return RESULT;
}
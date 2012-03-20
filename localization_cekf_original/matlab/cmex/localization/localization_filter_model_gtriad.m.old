%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var_out = localization_filter_model_gtriad(procedure_name, X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias, flagconsiderprediction)

q_predicted = X(1:4,1);
switch procedure_name
    case 'evaluate'
        if flagestimateaccelerometerbias
            bax = X(11); bay = X(12); baz = X(13);
        else
            bax = 0; bay = 0; baz = 0;
        end
        imumeasure.ax = imumeasure.ax - bax; imumeasure.ay = imumeasure.ay - bay; imumeasure.az = imumeasure.az - baz;
        if flagconsiderprediction
            var_out = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
        else
            var_out = localization_triad(imumeasure, magnetometermeasure, M, G);
        end
    case 'dg_du_imu'
        delta = 1e-5;
        if flagestimateaccelerometerbias
            bax = X(11); bay = X(12); baz = X(13);
        else
            bax = 0; bay = 0; baz = 0;
        end
        imumeasure.ax = imumeasure.ax - bax; imumeasure.ay = imumeasure.ay - bay; imumeasure.az = imumeasure.az - baz;

        if flagconsiderprediction
            y = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            imumeasure.ax = imumeasure.ax + delta;
            y_ax = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            imumeasure.ax = imumeasure.ax - delta;
            imumeasure.ay = imumeasure.ay + delta;
            y_ay = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            imumeasure.ay = imumeasure.ay - delta;
            imumeasure.az = imumeasure.az + delta;
            y_az = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            imumeasure.az = imumeasure.az - delta;
        else
            y = localization_triad(imumeasure, magnetometermeasure, M, G);
            imumeasure.ax = imumeasure.ax + delta;
            y_ax = localization_triad(imumeasure, magnetometermeasure, M, G);
            imumeasure.ax = imumeasure.ax - delta;
            imumeasure.ay = imumeasure.ay + delta;
            y_ay = localization_triad(imumeasure, magnetometermeasure, M, G);
            imumeasure.ay = imumeasure.ay - delta;
            imumeasure.az = imumeasure.az + delta;
            y_az = localization_triad(imumeasure, magnetometermeasure, M, G);
            imumeasure.az = imumeasure.az - delta;
        end
        
        dg_du_imu = [(y_ax-y)/delta (y_ay-y)/delta (y_az-y)/delta zeros(4,3)];
        var_out = dg_du_imu;
    case 'dg_du_mag'
        delta = 1e-5;
        if flagestimateaccelerometerbias
            bax = X(11); bay = X(12); baz = X(13);
        else
            bax = 0; bay = 0; baz = 0;
        end
        imumeasure.ax = imumeasure.ax - bax; imumeasure.ay = imumeasure.ay - bay; imumeasure.az = imumeasure.az - baz;

        if flagconsiderprediction
            y = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            magnetometermeasure.mx = magnetometermeasure.mx + delta;
            y_mx = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            magnetometermeasure.mx = magnetometermeasure.mx - delta;
            magnetometermeasure.my = magnetometermeasure.my + delta;
            y_my = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            magnetometermeasure.my = magnetometermeasure.my - delta;
            magnetometermeasure.mz = magnetometermeasure.mz + delta;
            y_mz = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
            magnetometermeasure.mz = magnetometermeasure.mz - delta;
        else
            y = localization_triad(imumeasure, magnetometermeasure, M, G);
            magnetometermeasure.mx = magnetometermeasure.mx + delta;
            y_mx = localization_triad(imumeasure, magnetometermeasure, M, G);
            magnetometermeasure.mx = magnetometermeasure.mx - delta;
            magnetometermeasure.my = magnetometermeasure.my + delta;
            y_my = localization_triad(imumeasure, magnetometermeasure, M, G);
            magnetometermeasure.my = magnetometermeasure.my - delta;
            magnetometermeasure.mz = magnetometermeasure.mz + delta;
            y_mz = localization_triad(imumeasure, magnetometermeasure, M, G);
            magnetometermeasure.mz = magnetometermeasure.mz - delta;
        end
        dg_du_mag = [(y_mx-y)/delta (y_my-y)/delta (y_mz-y)/delta];
        var_out = dg_du_mag;

    case 'dg_dx'
        delta = 1e-5;
        if flagestimateaccelerometerbias
            bax = X(11); bay = X(12); baz = X(13);
            imumeasure.ax = imumeasure.ax - bax; imumeasure.ay = imumeasure.ay - bay; imumeasure.az = imumeasure.az - baz;

            if flagconsiderprediction
                y = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
                imumeasure.ax = imumeasure.ax - delta;
                y_ax = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
                imumeasure.ax = imumeasure.ax + delta;
                imumeasure.ay = imumeasure.ay - delta;
                y_ay = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
                imumeasure.ay = imumeasure.ay + delta;
                imumeasure.az = imumeasure.az - delta;
                y_az = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
                imumeasure.az = imumeasure.az + delta;
            else
                y = localization_triad(imumeasure, magnetometermeasure, M, G);
                imumeasure.ax = imumeasure.ax - delta;
                y_ax = localization_triad(imumeasure, magnetometermeasure, M, G);
                imumeasure.ax = imumeasure.ax + delta;
                imumeasure.ay = imumeasure.ay - delta;
                y_ay = localization_triad(imumeasure, magnetometermeasure, M, G);
                imumeasure.ay = imumeasure.ay + delta;
                imumeasure.az = imumeasure.az - delta;
                y_az = localization_triad(imumeasure, magnetometermeasure, M, G);
                imumeasure.az = imumeasure.az + delta;
            end
            dg_dx = [zeros(4,10) (y_ax-y)/delta (y_ay-y)/delta (y_az-y)/delta];
        else
            dg_dx = zeros(4,length(X));
        end
        var_out = dg_dx;
end
<img src="https://github.com/mitsuba-renderer/mitsuba2/raw/master/docs/images/logo_plain.png" width="120" height="120" alt="Mitsuba logo">

# Mitsuba 2-NLOS

**Mitsuba 2, extended for transient path tracing and non-line of sight data capture.**

## Usage

* Transient path tracing: see [transientintegrator](https://github.com/diegoroyo/mitsuba2/blob/feat-transient/src/librender/transientintegrator.cpp), [transientpath](https://github.com/diegoroyo/mitsuba2/blob/feat-transient/src/integrators/transientpath.cpp) and [transientstokes](https://github.com/diegoroyo/mitsuba2/blob/feat-transient/src/integrators/transientstokes.cpp).
* Non-line of sight data capture: see [nloscapturemeter](https://github.com/diegoroyo/mitsuba2/blob/feat-transient/src/sensors/nloscapturemeter.cpp) (see documentation)

<details>
   <summary>Example scene for line-of-sight capture</summary>

   ```xml
   <scene version="2.2.1">
      <!-- Use transientpath integrator -->
      <integrator type="transientpath">
         <!-- Discard paths with depth >= max_depth -->
         <integer name="max_depth" value="4"/>
      </integrator>

      <!-- Geometry, sensor, etc. -->
   </scene>
   ```
</details>

<details>
   <summary>Example scene for non-line of sight capture</summary>

   _Note that variables that start with $ should be changed_

   ```xml
   <scene version="2.2.1">
      <integrator type="transientpath">
         <!-- Recommended 1 for progress bar (see path integrator) -->
         <integer name="block_size" value="1"/>
         <!-- Discard paths with depth >= max_depth -->
         <integer name="max_depth" value="4"/>
         <!-- Only account for paths with depth = filter_depth -->
         <!-- <integer name="filter_depth" value="3"/> -->
         <boolean name="discard_direct_paths" value="true"/>
         <!-- Next event estimation for the laser through the relay wall (recommended true) -->
         <boolean name="nlos_emitter_sampling" value="true"/>
      </integrator>

      <!-- Relay wall and hidden geometry materials -->
      <bsdf type="diffuse" id="white">
         <rgb name="reflectance" value="0.7, 0.7, 0.7"/>
      </bsdf>

      <!-- Hidden geometry -->
      <shape type="obj">
         <string name="filename" value="$hidden_mesh_obj"/>
         <bsdf type="twosided">
               <ref id="white"/>
         </bsdf>

         <transform name="to_world">
               <scale x="$hidden_scale" y="$hidden_scale" z="$hidden_scale"/>
               <rotate x="1" angle="$hidden_rot_degrees_x"/>
               <rotate y="1" angle="$hidden_rot_degrees_y"/>
               <rotate z="1" angle="$hidden_rot_degrees_z"/>
               <translate x="0" y="0" z="$hidden_distance_to_wall"/>
         </transform>
      </shape>

      <!-- Relay wall -->
      <shape type="rectangle">
         <ref id="white"/>

         <transform name="to_world">
               <rotate z="1" angle="180"/>
               <scale x="$relay_wall_scale" y="$relay_wall_scale" z="1"/>
         </transform>

         <!-- NLOS capture sensor placed on the relay wall -->
         <sensor type="nloscapturemeter">
               <sampler type="independent">
                  <integer name="sample_count" value="$sample_count"/>
               </sampler>

               <!-- Laser -->
               <emitter type="projector">
                  <spectrum name="irradiance" value="400:0, 500:80, 600:156.0, 700:184.0"/>
                  <float name="fov" value="1.5"/>
               </emitter>

               <!-- Acount time of flight for the laser->relay wall and relay wall->sensor paths -->
               <boolean name="account_first_and_last_bounces" value="$account_first_and_last_bounces"/>

               <!-- World-space coordinates -->
               <point name="sensor_origin" x="-0.5" y="0" z="0.25"/>
               <point name="laser_origin" x="-0.5" y="0" z="0.25"/>
               
               <!-- alternative to laser_lookat_pixel -->
               <!-- <point name="laser_lookat_3d" x="0" y="0" z="0"/> -->
               
               <!-- Screen-space coordinates (see streakhdrfilm) -->
               <point name="laser_lookat_pixel" x="$laser_lookat_x" y="$laser_lookat_y" z="0"/>

               <!-- Transient image I(width, height, num_bins) -->
               <film type="streakhdrfilm" name="streakfilm">
                  <integer name="width" value="$sensor_width"/>
                  <integer name="height" value="$sensor_height"/>

                  <!-- Recommended to prevent clamping -->
                  <string name="component_format" value="float32"/>

                  <integer name="num_bins" value="$num_bins"/>

                  <!-- Auto-detect start_opl (and also bin_width_opl if set to a negative value) -->
                  <boolean name="auto_detect_bins" value="$auto_detect_bins"/>
                  <float name="bin_width_opl" value="$bin_width_opl"/>
                  <float name="start_opl" value="$start_opl"/>

                  <rfilter name="rfilter" type="box"/>
                  <!-- NOTE: tfilters are not implemented yet -->
                  <!-- <rfilter name="tfilter" type="box"/>  -->
                  <boolean name="high_quality_edges" value="false"/>
               </film>
         </sensor>
      </shape>
   </scene>
   ```
</details>

<br>

Refer to the [original Mitsuba 2 repository](https://github.com/mitsuba-renderer/mitsuba2)
for additional documentation and instructions on how to compile, use, and extend Mitsuba 2.

## License

See the [original repository](https://github.com/mitsuba-renderer/mitsuba2).
<!-- Additionally, if you are using this code in academic research, we would be grateful
if you cited our paper: -->

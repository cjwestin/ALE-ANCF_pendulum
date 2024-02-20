function Q = ALEANCF_axialforce(axial_stiffness,p1,p2,qe)

% Function to calculate generalized axial force 
% for an ALE-ANCF element
%
% Auto-generated using the Matlab Symbolic Toolbox
%
% INPUTS:
% axial_stiffness - cable axial stiffness EA (N)
% p1 - material coordinate of first element node
% p2 - material coordinate of second element node
% qe - vector of reduced generalized coordinates
%
% OUTPUTS:
% J - Jacobian matrix
%

qe1 = qe(1);
qe2 = qe(2);
qe3 = qe(3);
qe4 = qe(4);
qe5 = qe(5);
qe6 = qe(6);

qe7 = qe(7);
qe8 = qe(8);
qe9 = qe(9);
qe10 = qe(10);
qe11 = qe(11);
qe12 = qe(12);


t2 = qe1.*qe4;
t3 = qe2.*qe5;
t4 = qe3.*qe6;
t5 = qe1.*qe10;
t6 = qe4.*qe7;
t7 = qe2.*qe11;
t8 = qe5.*qe8;
t9 = qe3.*qe12;
t10 = qe6.*qe9;
t11 = qe7.*qe10;
t12 = qe8.*qe11;
t13 = qe9.*qe12;
t14 = qe1.^2;
t15 = qe2.^2;
t16 = qe3.^2;
t17 = qe4.^2;
t18 = qe4.^3;
t19 = qe5.^2;
t21 = qe5.^3;
t22 = qe6.^2;
t24 = qe6.^3;
t25 = qe7.^2;
t27 = qe8.^2;
t28 = qe9.^2;
t29 = qe10.^2;
t30 = qe10.^3;
t31 = qe11.^2;
t33 = qe11.^3;
t34 = qe12.^2;
t36 = qe12.^3;
t38 = qe1.*qe7.*2.0;
t39 = qe2.*qe8.*2.0;
t40 = qe3.*qe9.*2.0;
t41 = -p2;
t42 = -qe7;
t43 = -qe8;
t44 = -qe9;
t45 = qe1.*1.4e+1;
t46 = qe2.*1.4e+1;
t47 = qe3.*1.4e+1;
t48 = qe4.*1.4e+1;
t49 = qe5.*1.4e+1;
t50 = qe6.*1.4e+1;
t51 = qe7.*1.4e+1;
t52 = qe8.*1.4e+1;
t53 = qe9.*1.4e+1;
t54 = qe10.*1.4e+1;
t55 = qe11.*1.4e+1;
t56 = qe12.*1.4e+1;
t57 = qe1.*4.2e+1;
t58 = qe2.*4.2e+1;
t59 = qe3.*4.2e+1;
t60 = qe7.*4.2e+1;
t61 = qe8.*4.2e+1;
t62 = qe9.*4.2e+1;
t75 = qe1.*qe7.*2.8e+1;
t76 = qe2.*qe8.*2.8e+1;
t77 = qe3.*qe9.*2.8e+1;
t200 = qe4.*qe5.*qe10.*2.0;
t209 = qe2.*qe4.*qe10.*6.0;
t214 = qe4.*qe5.*qe11.*2.0;
t215 = qe4.*qe6.*qe10.*2.0;
t220 = qe1.*qe5.*qe11.*6.0;
t224 = qe3.*qe4.*qe10.*6.0;
t236 = qe4.*qe6.*qe12.*2.0;
t237 = qe5.*qe6.*qe11.*2.0;
t240 = qe1.*qe6.*qe12.*6.0;
t247 = qe3.*qe5.*qe11.*6.0;
t253 = qe5.*qe6.*qe12.*2.0;
t255 = qe2.*qe6.*qe12.*6.0;
t262 = qe4.*qe10.*qe11.*2.0;
t267 = qe4.*qe8.*qe10.*6.0;
t268 = qe4.*qe10.*qe12.*2.0;
t270 = qe5.*qe10.*qe11.*2.0;
t277 = qe4.*qe9.*qe10.*6.0;
t278 = qe5.*qe7.*qe11.*6.0;
t282 = qe5.*qe11.*qe12.*2.0;
t283 = qe6.*qe10.*qe12.*2.0;
t291 = qe5.*qe9.*qe11.*6.0;
t292 = qe6.*qe7.*qe12.*6.0;
t295 = qe6.*qe11.*qe12.*2.0;
t301 = qe6.*qe8.*qe12.*6.0;
t341 = qe4.*qe5.*qe10.*qe11.*4.0;
t344 = qe4.*qe6.*qe10.*qe12.*4.0;
t347 = qe5.*qe6.*qe11.*qe12.*4.0;
t363 = qe4.*qe10.*-1.4e+1;
t364 = qe5.*qe11.*-1.4e+1;
t365 = qe6.*qe12.*-1.4e+1;
t434 = qe1.*qe5.*qe7.*-2.0;
t440 = qe1.*qe6.*qe7.*-2.0;
t441 = qe2.*qe4.*qe8.*-2.0;
t448 = qe2.*qe6.*qe8.*-2.0;
t449 = qe3.*qe4.*qe9.*-2.0;
t453 = qe3.*qe5.*qe9.*-2.0;
t455 = qe1.*qe7.*qe11.*-2.0;
t461 = qe1.*qe7.*qe12.*-2.0;
t466 = qe2.*qe8.*qe10.*-2.0;
t477 = qe2.*qe8.*qe12.*-2.0;
t482 = qe3.*qe9.*qe10.*-2.0;
t493 = qe3.*qe9.*qe11.*-2.0;
t20 = t17.^2;
t23 = t19.^2;
t26 = t22.^2;
t32 = t29.^2;
t35 = t31.^2;
t37 = t34.^2;
t63 = -t38;
t64 = -t39;
t65 = -t6;
t66 = -t40;
t67 = -t8;
t68 = -t10;
t69 = -t11;
t70 = -t12;
t71 = -t13;
t72 = qe10.*t48;
t73 = qe11.*t49;
t74 = qe12.*t50;
t78 = qe1.*t19;
t79 = qe2.*t17;
t80 = qe5.*t14;
t81 = qe4.*t15;
t82 = qe1.*t22;
t83 = qe3.*t17;
t84 = qe6.*t14;
t85 = qe4.*t16;
t86 = qe2.*t22;
t87 = qe3.*t19;
t88 = qe6.*t15;
t89 = qe5.*t16;
t90 = qe4.*t19;
t91 = qe5.*t17;
t92 = qe4.*t22;
t93 = qe6.*t17;
t94 = qe5.*t22;
t95 = qe6.*t19;
t96 = qe1.*t31;
t97 = qe2.*t29;
t98 = qe4.*t27;
t99 = qe5.*t25;
t100 = qe11.*t14;
t101 = qe10.*t15;
t102 = qe8.*t17;
t103 = qe7.*t19;
t104 = qe1.*t34;
t105 = qe3.*t29;
t106 = qe4.*t28;
t107 = qe6.*t25;
t108 = qe12.*t14;
t109 = qe10.*t16;
t110 = qe9.*t17;
t111 = qe7.*t22;
t112 = qe2.*t34;
t113 = qe3.*t31;
t114 = qe5.*t28;
t115 = qe6.*t27;
t116 = qe12.*t15;
t117 = qe11.*t16;
t118 = qe9.*t19;
t119 = qe8.*t22;
t120 = qe4.*t31;
t121 = qe5.*t29;
t122 = qe11.*t17;
t123 = qe10.*t19;
t124 = qe4.*t34;
t125 = qe6.*t29;
t126 = qe12.*t17;
t127 = qe10.*t22;
t128 = qe5.*t34;
t129 = qe6.*t31;
t130 = qe12.*t19;
t131 = qe11.*t22;
t132 = qe7.*t31;
t133 = qe8.*t29;
t134 = qe11.*t25;
t135 = qe10.*t27;
t136 = qe7.*t34;
t137 = qe9.*t29;
t138 = qe12.*t25;
t139 = qe10.*t28;
t140 = qe8.*t34;
t141 = qe9.*t31;
t142 = qe12.*t27;
t143 = qe11.*t28;
t144 = qe10.*t31;
t145 = qe11.*t29;
t146 = qe10.*t34;
t147 = qe12.*t29;
t148 = qe11.*t34;
t149 = qe12.*t31;
t150 = qe2.*t2.*2.0;
t151 = qe1.*t3.*2.0;
t152 = qe3.*t2.*2.0;
t153 = qe1.*t4.*2.0;
t154 = qe5.*t2.*2.0;
t155 = qe3.*t3.*2.0;
t156 = qe6.*t2.*2.0;
t157 = qe2.*t4.*2.0;
t158 = qe4.*t3.*2.0;
t159 = qe2.*t5.*2.0;
t160 = qe8.*t2.*2.0;
t161 = qe5.*t38;
t162 = qe2.*t6.*2.0;
t163 = qe6.*t3.*2.0;
t164 = qe4.*t4.*2.0;
t165 = qe1.*t7.*2.0;
t166 = qe3.*t5.*2.0;
t167 = qe5.*t2.*6.0;
t168 = qe9.*t2.*2.0;
t169 = qe1.*t8.*2.0;
t170 = qe6.*t38;
t171 = qe4.*t39;
t172 = qe7.*t3.*2.0;
t173 = qe3.*t6.*2.0;
t174 = qe5.*t4.*2.0;
t175 = qe6.*t2.*6.0;
t176 = qe4.*t3.*6.0;
t177 = qe1.*t9.*2.0;
t178 = qe7.*t2.*6.0;
t179 = qe1.*t10.*2.0;
t180 = qe3.*t7.*2.0;
t181 = qe9.*t3.*2.0;
t182 = qe6.*t39;
t183 = qe4.*t40;
t184 = qe3.*t8.*2.0;
t185 = qe7.*t4.*2.0;
t186 = qe5.*t6.*2.0;
t187 = qe2.*t9.*2.0;
t188 = qe6.*t3.*6.0;
t189 = qe2.*t10.*2.0;
t190 = qe4.*t4.*6.0;
t191 = qe5.*t40;
t192 = qe8.*t4.*2.0;
t193 = qe4.*t8.*2.0;
t194 = qe6.*t6.*2.0;
t195 = qe5.*t4.*6.0;
t196 = qe11.*t38;
t197 = qe8.*t5.*2.0;
t198 = qe8.*t3.*6.0;
t199 = qe2.*t11.*2.0;
t201 = qe4.*t10.*2.0;
t202 = qe8.*t6.*2.0;
t203 = qe6.*t8.*2.0;
t204 = qe11.*t2.*6.0;
t205 = qe5.*t5.*6.0;
t206 = qe12.*t38;
t207 = qe1.*t12.*2.0;
t208 = qe9.*t5.*2.0;
t210 = qe7.*t7.*2.0;
t211 = qe10.*t39;
t212 = qe3.*t11.*2.0;
t213 = qe5.*t6.*6.0;
t216 = qe9.*t6.*2.0;
t217 = qe5.*t10.*2.0;
t218 = qe7.*t8.*2.0;
t219 = qe12.*t2.*6.0;
t221 = qe6.*t5.*6.0;
t222 = qe4.*t7.*6.0;
t223 = qe10.*t3.*6.0;
t225 = qe4.*t8.*6.0;
t226 = qe6.*t6.*6.0;
t227 = qe7.*t5.*6.0;
t228 = qe1.*t13.*2.0;
t229 = qe11.*t5.*2.0;
t230 = qe12.*t39;
t231 = qe9.*t7.*2.0;
t232 = qe9.*t4.*6.0;
t233 = qe7.*t9.*2.0;
t234 = qe3.*t12.*2.0;
t235 = qe10.*t40;
t238 = qe9.*t8.*2.0;
t239 = qe7.*t10.*2.0;
t241 = qe12.*t5.*2.0;
t242 = qe12.*t3.*6.0;
t243 = qe6.*t7.*6.0;
t244 = qe2.*t13.*2.0;
t245 = qe10.*t7.*2.0;
t246 = qe4.*t9.*6.0;
t248 = qe10.*t4.*6.0;
t249 = qe8.*t9.*2.0;
t250 = qe11.*t40;
t251 = qe4.*t10.*6.0;
t252 = qe6.*t8.*6.0;
t254 = qe8.*t10.*2.0;
t256 = qe5.*t9.*6.0;
t257 = qe11.*t4.*6.0;
t258 = qe5.*t10.*6.0;
t259 = qe8.*t7.*6.0;
t260 = qe12.*t7.*2.0;
t261 = qe10.*t9.*2.0;
t263 = qe8.*t11.*2.0;
t264 = qe11.*t5.*6.0;
t265 = qe11.*t9.*2.0;
t266 = qe11.*t6.*6.0;
t269 = qe5.*t11.*6.0;
t271 = qe7.*t12.*2.0;
t272 = qe9.*t11.*2.0;
t273 = qe12.*t5.*6.0;
t274 = qe10.*t7.*6.0;
t275 = qe12.*t6.*6.0;
t276 = qe4.*t12.*6.0;
t279 = qe10.*t8.*6.0;
t280 = qe6.*t11.*6.0;
t281 = qe9.*t9.*6.0;
t284 = qe7.*t13.*2.0;
t285 = qe11.*t11.*2.0;
t286 = qe9.*t12.*2.0;
t287 = qe12.*t7.*6.0;
t288 = qe10.*t9.*6.0;
t289 = qe4.*t13.*6.0;
t290 = qe12.*t8.*6.0;
t293 = qe6.*t12.*6.0;
t294 = qe10.*t10.*6.0;
t296 = qe12.*t11.*2.0;
t297 = qe8.*t13.*2.0;
t298 = qe10.*t12.*2.0;
t299 = qe11.*t9.*6.0;
t300 = qe5.*t13.*6.0;
t302 = qe11.*t10.*6.0;
t303 = qe12.*t12.*2.0;
t304 = qe10.*t13.*2.0;
t305 = qe11.*t11.*6.0;
t306 = qe11.*t13.*2.0;
t307 = qe12.*t11.*6.0;
t308 = qe10.*t12.*6.0;
t309 = qe12.*t12.*6.0;
t310 = qe10.*t13.*6.0;
t311 = qe11.*t13.*6.0;
t312 = -t48;
t313 = -t49;
t314 = -t50;
t315 = -t51;
t316 = -t52;
t317 = -t53;
t318 = -t54;
t319 = -t55;
t320 = -t56;
t321 = -t60;
t322 = -t61;
t323 = -t62;
t324 = t2.*t3.*4.0;
t325 = t2.*t4.*4.0;
t326 = t3.*t4.*4.0;
t327 = t2.*t8.*4.0;
t328 = t3.*t6.*4.0;
t329 = t2.*t10.*4.0;
t330 = t4.*t6.*4.0;
t331 = t3.*t10.*4.0;
t332 = t4.*t8.*4.0;
t333 = t5.*t7.*4.0;
t334 = t6.*t8.*4.0;
t335 = t5.*t9.*4.0;
t336 = t6.*t10.*4.0;
t337 = t7.*t9.*4.0;
t338 = t8.*t10.*4.0;
t339 = t5.*t12.*4.0;
t340 = t7.*t11.*4.0;
t342 = t5.*t13.*4.0;
t343 = t9.*t11.*4.0;
t345 = t7.*t13.*4.0;
t346 = t9.*t12.*4.0;
t348 = t11.*t12.*4.0;
t349 = t11.*t13.*4.0;
t350 = t12.*t13.*4.0;
t351 = -t18;
t352 = -t21;
t353 = -t24;
t354 = -t30;
t355 = -t33;
t356 = -t36;
t357 = t14.*1.4e+1;
t358 = t15.*1.4e+1;
t359 = t16.*1.4e+1;
t360 = t25.*1.4e+1;
t361 = t27.*1.4e+1;
t362 = t28.*1.4e+1;
t366 = p1+t41;
t367 = qe1+t42;
t368 = qe2+t43;
t369 = qe3+t44;
t370 = qe4.*t2.*3.0;
t371 = qe1.*t2.*3.0;
t375 = qe5.*t3.*3.0;
t377 = qe2.*t3.*3.0;
t380 = qe6.*t4.*3.0;
t381 = qe3.*t4.*3.0;
t382 = qe10.*t5.*3.0;
t383 = qe7.*t6.*3.0;
t384 = qe1.*t5.*3.0;
t385 = qe4.*t6.*3.0;
t391 = qe11.*t7.*3.0;
t393 = qe8.*t8.*3.0;
t394 = qe2.*t7.*3.0;
t396 = qe5.*t8.*3.0;
t400 = qe4.*t29.*3.0;
t401 = qe10.*t17.*3.0;
t404 = qe12.*t9.*3.0;
t405 = qe4.*t30.*3.0;
t406 = qe9.*t10.*3.0;
t407 = qe3.*t9.*3.0;
t408 = qe10.*t18.*3.0;
t409 = qe6.*t10.*3.0;
t410 = qe5.*t31.*3.0;
t411 = qe11.*t19.*3.0;
t412 = qe5.*t33.*3.0;
t413 = qe10.*t11.*3.0;
t414 = qe11.*t21.*3.0;
t415 = qe7.*t11.*3.0;
t416 = qe6.*t34.*3.0;
t419 = qe12.*t22.*3.0;
t420 = qe6.*t36.*3.0;
t422 = qe11.*t12.*3.0;
t424 = qe12.*t24.*3.0;
t425 = qe8.*t12.*3.0;
t428 = qe12.*t13.*3.0;
t429 = qe9.*t13.*3.0;
t464 = -t209;
t469 = -t220;
t473 = -t224;
t483 = -t240;
t490 = -t247;
t494 = -t255;
t502 = qe10.*t2.*1.8e+1;
t503 = qe11.*t3.*1.8e+1;
t504 = qe12.*t4.*1.8e+1;
t505 = qe10.*t6.*1.8e+1;
t506 = qe11.*t8.*1.8e+1;
t507 = qe12.*t10.*1.8e+1;
t508 = t14.*t19;
t509 = t15.*t17;
t510 = t14.*t22;
t511 = t16.*t17;
t512 = t15.*t22;
t513 = t16.*t19;
t514 = t14.*t31;
t515 = t15.*t29;
t516 = t17.*t27;
t517 = t19.*t25;
t518 = t14.*t34;
t519 = t16.*t29;
t520 = t17.*t28;
t521 = t22.*t25;
t522 = t15.*t34;
t523 = t16.*t31;
t524 = t19.*t28;
t525 = t22.*t27;
t526 = t17.*t31;
t527 = t19.*t29;
t528 = t17.*t34;
t529 = t22.*t29;
t530 = t19.*t34;
t531 = t22.*t31;
t532 = t25.*t31;
t533 = t27.*t29;
t534 = t25.*t34;
t535 = t28.*t29;
t536 = t27.*t34;
t537 = t28.*t31;
t538 = t19.*t38;
t539 = t22.*t38;
t540 = t17.*t39;
t541 = t2.*t6.*6.0;
t542 = t22.*t39;
t543 = t17.*t40;
t544 = t19.*t40;
t545 = t31.*t38;
t546 = t3.*t8.*6.0;
t547 = t34.*t38;
t548 = t29.*t39;
t552 = t5.*t11.*6.0;
t553 = t34.*t39;
t554 = t29.*t40;
t555 = t4.*t10.*6.0;
t556 = t31.*t40;
t560 = t7.*t12.*6.0;
t564 = t9.*t13.*6.0;
t576 = -t341;
t579 = -t344;
t582 = -t347;
t606 = qe4.*t2.*9.0;
t612 = qe5.*t3.*9.0;
t614 = qe6.*t4.*9.0;
t616 = qe10.*t5.*9.0;
t617 = qe4.*t6.*9.0;
t618 = qe11.*t7.*9.0;
t620 = qe5.*t8.*9.0;
t624 = qe12.*t9.*9.0;
t626 = qe6.*t10.*9.0;
t628 = qe10.*t11.*9.0;
t629 = qe11.*t12.*9.0;
t630 = qe12.*t13.*9.0;
t634 = t2.^2.*3.0;
t635 = t3.^2.*3.0;
t636 = t4.^2.*3.0;
t637 = t5.^2.*3.0;
t638 = t6.^2.*3.0;
t639 = t7.^2.*3.0;
t640 = t8.^2.*3.0;
t641 = t17.*t29.*3.0;
t642 = t9.^2.*3.0;
t643 = t10.^2.*3.0;
t644 = t19.*t31.*3.0;
t645 = t11.^2.*3.0;
t646 = t22.*t34.*3.0;
t647 = t12.^2.*3.0;
t648 = t13.^2.*3.0;
t372 = t78.*3.0;
t373 = t79.*3.0;
t374 = t82.*3.0;
t376 = t83.*3.0;
t378 = t86.*3.0;
t379 = t87.*3.0;
t386 = t96.*3.0;
t387 = t97.*3.0;
t388 = t102.*3.0;
t389 = t103.*3.0;
t390 = t104.*3.0;
t392 = t105.*3.0;
t395 = t110.*3.0;
t397 = t111.*3.0;
t398 = t112.*3.0;
t399 = t113.*3.0;
t402 = t118.*3.0;
t403 = t119.*3.0;
t417 = t132.*3.0;
t418 = t133.*3.0;
t421 = t136.*3.0;
t423 = t137.*3.0;
t426 = t140.*3.0;
t427 = t141.*3.0;
t430 = -t154;
t431 = -t156;
t432 = -t158;
t433 = -t160;
t435 = -t162;
t436 = -t163;
t437 = -t164;
t438 = -t168;
t439 = -t169;
t442 = -t172;
t443 = -t173;
t444 = -t174;
t445 = -t178;
t446 = -t179;
t447 = -t181;
t450 = -t184;
t451 = -t185;
t452 = -t189;
t454 = -t192;
t456 = -t197;
t457 = -t198;
t458 = -t199;
t459 = -t204;
t460 = -t205;
t462 = -t207;
t463 = -t208;
t465 = -t210;
t467 = -t212;
t468 = -t219;
t470 = -t221;
t471 = -t222;
t472 = -t223;
t474 = -t227;
t475 = -t228;
t476 = -t229;
t478 = -t231;
t479 = -t232;
t480 = -t233;
t481 = -t234;
t484 = -t241;
t485 = -t242;
t486 = -t243;
t487 = -t244;
t488 = -t245;
t489 = -t246;
t491 = -t248;
t492 = -t249;
t495 = -t256;
t496 = -t257;
t497 = -t259;
t498 = -t260;
t499 = -t261;
t500 = -t265;
t501 = -t281;
t549 = qe10.*t90.*3.0;
t550 = qe10.*t92.*3.0;
t551 = qe11.*t91.*3.0;
t557 = qe11.*t94.*3.0;
t558 = qe12.*t93.*3.0;
t559 = qe12.*t95.*3.0;
t561 = qe10.*t120.*3.0;
t562 = qe10.*t124.*3.0;
t563 = qe11.*t121.*3.0;
t565 = qe11.*t128.*3.0;
t566 = qe12.*t125.*3.0;
t567 = qe12.*t129.*3.0;
t568 = -t327;
t569 = -t328;
t570 = -t329;
t571 = -t330;
t572 = -t331;
t573 = -t332;
t574 = -t339;
t575 = -t340;
t577 = -t342;
t578 = -t343;
t580 = -t345;
t581 = -t346;
t583 = -t357;
t584 = -t358;
t585 = -t359;
t586 = -t360;
t587 = -t361;
t588 = -t362;
t589 = -t78;
t590 = -t79;
t591 = -t370;
t592 = -t82;
t593 = -t83;
t594 = -t86;
t595 = -t87;
t596 = -t375;
t597 = -t90;
t598 = -t91;
t599 = -t92;
t600 = -t93;
t601 = -t380;
t602 = -t94;
t603 = -t95;
t604 = -t96;
t605 = -t97;
t607 = -t382;
t608 = -t104;
t609 = -t105;
t610 = -t112;
t611 = -t113;
t613 = -t391;
t615 = -t404;
t619 = -t144;
t621 = -t145;
t622 = -t146;
t623 = -t147;
t625 = -t148;
t627 = -t149;
t631 = -t502;
t632 = -t503;
t633 = -t504;
t649 = qe7.*t78.*-2.0;
t650 = qe7.*t82.*-2.0;
t651 = qe8.*t79.*-2.0;
t652 = -t541;
t653 = qe8.*t86.*-2.0;
t654 = qe9.*t83.*-2.0;
t655 = qe9.*t87.*-2.0;
t656 = qe7.*t96.*-2.0;
t657 = -t546;
t658 = qe7.*t104.*-2.0;
t659 = qe8.*t97.*-2.0;
t660 = -t552;
t661 = qe8.*t112.*-2.0;
t662 = qe9.*t105.*-2.0;
t663 = -t555;
t664 = qe9.*t113.*-2.0;
t665 = -t560;
t666 = -t564;
t667 = -t526;
t668 = -t527;
t669 = -t641;
t670 = -t528;
t671 = -t529;
t672 = -t530;
t673 = -t531;
t674 = -t644;
t675 = -t646;
t676 = 1.0./t366;
t680 = t14+t15+t16+t25+t27+t28+t63+t64+t66;
t682 = t2+t3+t4+t5+t7+t9+t65+t67+t68+t69+t70+t71;
t677 = t676.^2;
t678 = t676.^3;
t681 = t680.^2;
t691 = t120+t123+t124+t127+t214+t236+t270+t283+t312+t318+t351+t354+t400+t401+t597+t599+t619+t622;
t692 = t121+t122+t128+t131+t200+t253+t262+t295+t313+t319+t352+t355+t410+t411+t598+t602+t621+t625;
t693 = t125+t126+t129+t130+t215+t237+t268+t282+t314+t320+t353+t356+t416+t419+t600+t603+t623+t627;
t697 = t45+t103+t111+t132+t136+t193+t201+t298+t304+t315+t385+t413+t432+t437+t488+t499+t589+t591+t592+t604+t607+t608;
t698 = t46+t102+t119+t133+t140+t186+t217+t285+t306+t316+t396+t422+t430+t444+t476+t500+t590+t594+t596+t605+t610+t613;
t699 = t47+t110+t118+t137+t141+t194+t203+t296+t303+t317+t409+t428+t431+t436+t484+t498+t593+t595+t601+t609+t611+t615;
t703 = t81+t85+t98+t101+t106+t109+t135+t139+t151+t153+t165+t177+t218+t239+t271+t284+t371+t383+t384+t415+t439+t441+t442+t445+t446+t449+t451+t462+t465+t466+t474+t475+t480+t482;
t704 = t80+t89+t99+t100+t114+t117+t134+t143+t150+t157+t159+t187+t202+t254+t263+t297+t377+t393+t394+t425+t433+t434+t435+t452+t453+t454+t455+t456+t457+t458+t487+t492+t493+t497;
t705 = t84+t88+t107+t108+t115+t116+t138+t142+t152+t155+t166+t180+t216+t238+t272+t286+t381+t406+t407+t429+t438+t440+t443+t447+t448+t450+t461+t463+t467+t477+t478+t479+t481+t501;
t709 = t75+t76+t77+t324+t325+t326+t333+t334+t335+t336+t337+t338+t348+t349+t350+t508+t509+t510+t511+t512+t513+t514+t515+t516+t517+t518+t519+t520+t521+t522+t523+t524+t525+t532+t533+t534+t535+t536+t537+t568+t569+t570+t571+t572+t573+t574+t575+t577+t578+t580+t581+t583+t584+t585+t586+t587+t588+t634+t635+t636+t637+t638+t639+t640+t642+t643+t645+t647+t648+t649+t650+t651+t652+t653+t654+t655+t656+t657+t658+t659+t660+t661+t662+t663+t664+t665+t666;
t679 = t677.^2;
t684 = axial_stiffness.*t367.*t678.*t680.*(3.6e+1./3.5e+1);
t685 = axial_stiffness.*t368.*t678.*t680.*(3.6e+1./3.5e+1);
t686 = axial_stiffness.*t369.*t678.*t680.*(3.6e+1./3.5e+1);
t687 = axial_stiffness.*t367.*t677.*t680.*(9.0./7.0e+1);
t688 = axial_stiffness.*t368.*t677.*t680.*(9.0./7.0e+1);
t689 = axial_stiffness.*t369.*t677.*t680.*(9.0./7.0e+1);
t690 = axial_stiffness.*t678.*t680.*t682.*(9.0./3.5e+1);
t694 = (axial_stiffness.*t691)./2.8e+2;
t695 = (axial_stiffness.*t692)./2.8e+2;
t696 = (axial_stiffness.*t693)./2.8e+2;
t700 = axial_stiffness.*t676.*t697.*(3.0./7.0e+1);
t701 = axial_stiffness.*t676.*t698.*(3.0./7.0e+1);
t702 = axial_stiffness.*t676.*t699.*(3.0./7.0e+1);
t706 = axial_stiffness.*t677.*t703.*(9.0./7.0e+1);
t707 = axial_stiffness.*t677.*t704.*(9.0./7.0e+1);
t708 = axial_stiffness.*t677.*t705.*(9.0./7.0e+1);
t710 = axial_stiffness.*t677.*t709.*(3.0./1.4e+2);
t683 = axial_stiffness.*t679.*t681.*(2.7e+1./3.5e+1);
Q = [-t684+t694+t700+t706;
        -t685+t695+t701+t707;
        -t686+t696+t702+t708;
        t687-(axial_stiffness.*(t57+t176+t190-t225-t251-t274+t276+t278+t279-t288+t289+t292+t294+t308+t310+t321+t372+t374-t386-t389-t390-t397+t417+t421+t469+t471+t472+t483+t489+t491+t505+t606-t616-t617+t628+t631-p1.*qe4.*5.6e+1+p2.*qe4.*5.6e+1-p2.*qe10.*1.4e+1+p1.*t18.*2.4e+1-p2.*t18.*2.4e+1-p1.*t30.*3.0+p2.*t30.*3.0+p1.*t54+p1.*t90.*2.4e+1-p2.*t90.*2.4e+1+p1.*t92.*2.4e+1-p2.*t92.*2.4e+1+p1.*t120.*2.0-p2.*t120.*2.0-p1.*t123.*3.0+p1.*t124.*2.0+p2.*t123.*3.0-p2.*t124.*2.0-p1.*t127.*3.0+p2.*t127.*3.0-p1.*t144.*3.0+p2.*t144.*3.0-p1.*t146.*3.0+p2.*t146.*3.0-p1.*qe10.*t17.*9.0+p2.*qe10.*t17.*9.0+p1.*qe4.*t29.*6.0-p2.*qe4.*t29.*6.0-p1.*qe4.*qe5.*qe11.*6.0+p2.*qe4.*qe5.*qe11.*6.0-p1.*qe4.*qe6.*qe12.*6.0+p2.*qe4.*qe6.*qe12.*6.0+p1.*qe5.*qe10.*qe11.*4.0-p2.*qe5.*qe10.*qe11.*4.0+p1.*qe6.*qe10.*qe12.*4.0-p2.*qe6.*qe10.*qe12.*4.0))./8.4e+2-axial_stiffness.*t676.*(t81+t85+t98+t106+t151+t153+t218+t239+t371+t383+t439+t441+t442+t445+t446+t449+t451).*(3.0./7.0e+1);
        t688-(axial_stiffness.*(t58+t167+t195-t213-t258-t264+t266+t267+t269-t299+t300+t301+t302+t305+t311+t322+t373+t378-t387-t388-t398-t403+t418+t426+t459+t460+t464+t494+t495+t496+t506+t612-t618-t620+t629+t632-p1.*qe5.*5.6e+1+p2.*qe5.*5.6e+1-p2.*qe11.*1.4e+1+p1.*t21.*2.4e+1-p2.*t21.*2.4e+1-p1.*t33.*3.0+p2.*t33.*3.0+p1.*t55+p1.*t91.*2.4e+1-p2.*t91.*2.4e+1+p1.*t94.*2.4e+1-p2.*t94.*2.4e+1+p1.*t121.*2.0-p1.*t122.*3.0-p2.*t121.*2.0+p2.*t122.*3.0+p1.*t128.*2.0-p2.*t128.*2.0-p1.*t131.*3.0+p2.*t131.*3.0-p1.*t145.*3.0+p2.*t145.*3.0-p1.*t148.*3.0+p2.*t148.*3.0-p1.*qe11.*t19.*9.0+p2.*qe11.*t19.*9.0+p1.*qe5.*t31.*6.0-p2.*qe5.*t31.*6.0-p1.*qe4.*qe5.*qe10.*6.0+p2.*qe4.*qe5.*qe10.*6.0-p1.*qe5.*qe6.*qe12.*6.0+p2.*qe5.*qe6.*qe12.*6.0+p1.*qe4.*qe10.*qe11.*4.0-p2.*qe4.*qe10.*qe11.*4.0+p1.*qe6.*qe11.*qe12.*4.0-p2.*qe6.*qe11.*qe12.*4.0))./8.4e+2-axial_stiffness.*t676.*(t80+t89+t99+t114+t150+t157+t202+t254+t377+t393+t433+t434+t435+t452+t453+t454+t457).*(3.0./7.0e+1);
        t689-(axial_stiffness.*(t59+t175+t188-t226-t252-t273+t275+t277+t280-t287+t290+t291+t293+t307+t309+t323+t376+t379-t392-t395-t399-t402+t423+t427+t468+t470+t473+t485+t486+t490+t507+t614-t624-t626+t630+t633-p1.*qe6.*5.6e+1+p2.*qe6.*5.6e+1-p2.*qe12.*1.4e+1+p1.*t24.*2.4e+1-p2.*t24.*2.4e+1-p1.*t36.*3.0+p2.*t36.*3.0+p1.*t56+p1.*t93.*2.4e+1-p2.*t93.*2.4e+1+p1.*t95.*2.4e+1-p2.*t95.*2.4e+1+p1.*t125.*2.0-p1.*t126.*3.0-p2.*t125.*2.0+p2.*t126.*3.0+p1.*t129.*2.0-p1.*t130.*3.0-p2.*t129.*2.0+p2.*t130.*3.0-p1.*t147.*3.0+p2.*t147.*3.0-p1.*t149.*3.0+p2.*t149.*3.0-p1.*qe12.*t22.*9.0+p2.*qe12.*t22.*9.0+p1.*qe6.*t34.*6.0-p2.*qe6.*t34.*6.0-p1.*qe4.*qe6.*qe10.*6.0+p2.*qe4.*qe6.*qe10.*6.0-p1.*qe5.*qe6.*qe11.*6.0+p2.*qe5.*qe6.*qe11.*6.0+p1.*qe4.*qe10.*qe12.*4.0-p2.*qe4.*qe10.*qe12.*4.0+p1.*qe5.*qe11.*qe12.*4.0-p2.*qe5.*qe11.*qe12.*4.0))./8.4e+2-axial_stiffness.*t676.*(t84+t88+t107+t115+t152+t155+t216+t238+t381+t406+t438+t440+t443+t447+t448+t450+t479).*(3.0./7.0e+1);
        t683-t690+t710+(axial_stiffness.*(t17.*-1.82e+2-t19.*1.82e+2+t20.*9.9e+1-t22.*1.82e+2+t23.*9.9e+1+t26.*9.9e+1+t29.*2.8e+1+t31.*2.8e+1-t32.*6.0+t34.*2.8e+1-t35.*6.0-t37.*6.0+t363+t364+t365+t405+t408+t412+t414+t420+t424+t549+t550+t551+t557+t558+t559+t561+t562+t563+t565+t566+t567+t576+t579+t582+t667+t668+t669+t670+t671+t672+t673+t674+t675+t17.*t19.*1.98e+2+t17.*t22.*1.98e+2+t19.*t22.*1.98e+2-t29.*t31.*1.2e+1-t29.*t34.*1.2e+1-t31.*t34.*1.2e+1))./8.4e+2;
        t684-t694-t700-t706;
        t685-t695-t701-t707;
        t686-t696-t702-t708;
        t687-(axial_stiffness.*(t57-t176-t190+t225+t251+t274+t276+t278+t279+t288+t289+t292+t294-t308-t310+t321-t372-t374+t386+t389+t390+t397-t417-t421+t469+t471+t472+t483+t489+t491+t505-t606+t616+t617-t628+t631-p2.*qe4.*1.4e+1-p1.*qe10.*5.6e+1+p2.*qe10.*5.6e+1-p1.*t18.*3.0+p2.*t18.*3.0+p1.*t30.*2.4e+1-p2.*t30.*2.4e+1+p1.*t48-p1.*t90.*3.0+p2.*t90.*3.0-p1.*t92.*3.0+p2.*t92.*3.0-p1.*t120.*3.0+p2.*t120.*3.0+p1.*t123.*2.0-p1.*t124.*3.0-p2.*t123.*2.0+p2.*t124.*3.0+p1.*t127.*2.0-p2.*t127.*2.0+p1.*t144.*2.4e+1-p2.*t144.*2.4e+1+p1.*t146.*2.4e+1-p2.*t146.*2.4e+1+p1.*qe10.*t17.*6.0-p2.*qe10.*t17.*6.0-p1.*qe4.*t29.*9.0+p2.*qe4.*t29.*9.0+p1.*qe4.*qe5.*qe11.*4.0-p2.*qe4.*qe5.*qe11.*4.0+p1.*qe4.*qe6.*qe12.*4.0-p2.*qe4.*qe6.*qe12.*4.0-p1.*qe5.*qe10.*qe11.*6.0+p2.*qe5.*qe10.*qe11.*6.0-p1.*qe6.*qe10.*qe12.*6.0+p2.*qe6.*qe10.*qe12.*6.0))./8.4e+2-axial_stiffness.*t676.*(t101+t109+t135+t139+t165+t177+t271+t284+t384+t415+t462+t465+t466+t474+t475+t480+t482).*(3.0./7.0e+1);
        t688-(axial_stiffness.*(t58-t167-t195+t213+t258+t264+t266+t267+t269+t299+t300+t301+t302-t305-t311+t322-t373-t378+t387+t388+t398+t403-t418-t426+t459+t460+t464+t494+t495+t496+t506-t612+t618+t620-t629+t632-p2.*qe5.*1.4e+1-p1.*qe11.*5.6e+1+p2.*qe11.*5.6e+1-p1.*t21.*3.0+p2.*t21.*3.0+p1.*t33.*2.4e+1-p2.*t33.*2.4e+1+p1.*t49-p1.*t91.*3.0+p2.*t91.*3.0-p1.*t94.*3.0+p2.*t94.*3.0-p1.*t121.*3.0+p1.*t122.*2.0+p2.*t121.*3.0-p2.*t122.*2.0-p1.*t128.*3.0+p2.*t128.*3.0+p1.*t131.*2.0-p2.*t131.*2.0+p1.*t145.*2.4e+1-p2.*t145.*2.4e+1+p1.*t148.*2.4e+1-p2.*t148.*2.4e+1+p1.*qe11.*t19.*6.0-p2.*qe11.*t19.*6.0-p1.*qe5.*t31.*9.0+p2.*qe5.*t31.*9.0+p1.*qe4.*qe5.*qe10.*4.0-p2.*qe4.*qe5.*qe10.*4.0+p1.*qe5.*qe6.*qe12.*4.0-p2.*qe5.*qe6.*qe12.*4.0-p1.*qe4.*qe10.*qe11.*6.0+p2.*qe4.*qe10.*qe11.*6.0-p1.*qe6.*qe11.*qe12.*6.0+p2.*qe6.*qe11.*qe12.*6.0))./8.4e+2-axial_stiffness.*t676.*(t100+t117+t134+t143+t159+t187+t263+t297+t394+t425+t455+t456+t458+t487+t492+t493+t497).*(3.0./7.0e+1);
        t689-(axial_stiffness.*(t59-t175-t188+t226+t252+t273+t275+t277+t280+t287+t290+t291+t293-t307-t309+t323-t376-t379+t392+t395+t399+t402-t423-t427+t468+t470+t473+t485+t486+t490+t507-t614+t624+t626-t630+t633-p2.*qe6.*1.4e+1-p1.*qe12.*5.6e+1+p2.*qe12.*5.6e+1-p1.*t24.*3.0+p2.*t24.*3.0+p1.*t36.*2.4e+1-p2.*t36.*2.4e+1+p1.*t50-p1.*t93.*3.0+p2.*t93.*3.0-p1.*t95.*3.0+p2.*t95.*3.0-p1.*t125.*3.0+p1.*t126.*2.0+p2.*t125.*3.0-p2.*t126.*2.0-p1.*t129.*3.0+p1.*t130.*2.0+p2.*t129.*3.0-p2.*t130.*2.0+p1.*t147.*2.4e+1-p2.*t147.*2.4e+1+p1.*t149.*2.4e+1-p2.*t149.*2.4e+1+p1.*qe12.*t22.*6.0-p2.*qe12.*t22.*6.0-p1.*qe6.*t34.*9.0+p2.*qe6.*t34.*9.0+p1.*qe4.*qe6.*qe10.*4.0-p2.*qe4.*qe6.*qe10.*4.0+p1.*qe5.*qe6.*qe11.*4.0-p2.*qe5.*qe6.*qe11.*4.0-p1.*qe4.*qe10.*qe12.*6.0+p2.*qe4.*qe10.*qe12.*6.0-p1.*qe5.*qe11.*qe12.*6.0+p2.*qe5.*qe11.*qe12.*6.0))./8.4e+2-axial_stiffness.*t676.*(t108+t116+t138+t142+t166+t180+t272+t286+t407+t429+t461+t463+t467+t477+t478+t481+t501).*(3.0./7.0e+1);
        -t683+t690-t710-(axial_stiffness.*(t17.*2.8e+1+t19.*2.8e+1-t20.*6.0+t22.*2.8e+1-t23.*6.0-t26.*6.0-t29.*1.82e+2-t31.*1.82e+2+t32.*9.9e+1-t34.*1.82e+2+t35.*9.9e+1+t37.*9.9e+1+t363+t364+t365+t405+t408+t412+t414+t420+t424+t549+t550+t551+t557+t558+t559+t561+t562+t563+t565+t566+t567+t576+t579+t582+t667+t668+t669+t670+t671+t672+t673+t674+t675-t17.*t19.*1.2e+1-t17.*t22.*1.2e+1-t19.*t22.*1.2e+1+t29.*t31.*1.98e+2+t29.*t34.*1.98e+2+t31.*t34.*1.98e+2))./8.4e+2];
�	2��̃@2��̃@!2��̃@	����{V@����{V@!����{V@"e
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails$2��̃@٘��l�?A|}�K�P@Y�=�4a��?*	�S㥛�\@2t
=Iterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate{����?!���LWB@)�;�_�E�?17�b^D@@:Preprocessing2j
3Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat�DIH�m�?!rK����7@)��c> Й?1��ն�5@:Preprocessing2F
Iterator::Model'���K�?!�X�Z�@@)�;�2T�?1T���� 2@:Preprocessing2S
Iterator::Model::ParallelMapI�0e��?!>C�_l},@)I�0e��?1>C�_l},@:Preprocessing2X
!Iterator::Model::ParallelMap::Zip1zn�+�?!�S�R �P@)/��0x?1 eҲ?L@:Preprocessing2�
MIterator::Model::ParallelMap::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice'���s?!<��t�@)'���s?1<��t�@:Preprocessing2v
?Iterator::Model::ParallelMap::Zip[1]::ForeverRepeat::FromTensor~�*O �d?!�����@)~�*O �d?1�����@:Preprocessing2d
-Iterator::Model::ParallelMap::Zip[0]::FlatMapf����?!���9NWC@)��UJ��b?1�>z0  @:Preprocessing:�
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
�Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
�Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
�Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
�Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)�
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis�
device�Your program is NOT input-bound because only 2.3% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no#You may skip the rest of this page.B�
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown�
	٘��l�?٘��l�?!٘��l�?      ��!       "      ��!       *      ��!       2	|}�K�P@|}�K�P@!|}�K�P@:      ��!       B      ��!       J	�=�4a��?�=�4a��?!�=�4a��?R      ��!       Z	�=�4a��?�=�4a��?!�=�4a��?JCPU_ONLY2black"�
device�Your program is NOT input-bound because only 2.3% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.b
`input_pipeline_analyzer (especially Section 3 for the breakdown of input operations on the Host)m
ktrace_viewer (look at the activities on the timeline of each Host Thread near the bottom of the trace view)"T
Rtensorflow_stats (identify the time-consuming operations executed on the CPU_ONLY)"Z
Xtrace_viewer (look at the activities on the timeline of each CPU_ONLY in the trace view)*y
w<a href="https://www.tensorflow.org/guide/data_performance" target="_blank">Better performance with the tf.data API</a>2I
=type.googleapis.com/tensorflow.profiler.GenericRecommendation
nono:
Refer to the TF2 Profiler FAQ2"CPU: 
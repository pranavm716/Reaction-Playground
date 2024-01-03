<template>
  <v-container class="ma-4 pa-4">
    <h1 :class="`d-inline-block text-h4 font-weight-medium text-${$route.meta.navBarColor}`">{{ highlightedWord }}</h1>
    <h1 class="d-inline-block text-h4 font-weight-light">{{ otherWords }}</h1>
  </v-container>
</template>

<script>

import { defineComponent } from 'vue'

export default defineComponent({
  name: 'PageHeader',
  props: ['highlightedWord', 'otherWords']
})
</script>
